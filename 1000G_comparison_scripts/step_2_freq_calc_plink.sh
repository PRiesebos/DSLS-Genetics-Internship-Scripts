#!/bin/bash

# Define variables
output_dir='/scratch/hb-functionalgenomics/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/SRP129004/filter_comparison'
output_vcf_dir="$output_dir/vcf"
freq_dir="$output_dir/freq"

# Step 4: Run plink command on output vcf files in parallel
for file in "$output_vcf_dir"/*.vcf.gz; do
    sbatch << EOF
#!/bin/bash
#SBATCH --job-name=step4_plink_${file##*/}
#SBATCH --output=slurm_logs/step4_plink_${file##*/}_output.txt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=16:00:00
#SBATCH --get-user-env=L

# Run plink command and output frequency files to freq directory
plink --vcf "$file" --freq --out "$freq_dir/$(basename "${file%.vcf.gz}")"
EOF
done

# End of the script

