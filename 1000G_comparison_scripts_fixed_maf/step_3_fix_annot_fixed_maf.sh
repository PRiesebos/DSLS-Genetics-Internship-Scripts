#!/bin/bash

# Define variables
output_dir='/scratch/hb-functionalgenomics/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/SRP129004/filter_comparison_fixed_maf'
freq_dir="$output_dir/freq"
fixed_freq_dir="$output_dir/fixed_freq"
fix_annot_script_dir='/scratch/hb-functionalgenomics/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/extra_scripts/fix_annot.py'

# Step 5: Run fix_annot script on output frq files in parallel
for frq_file in "$freq_dir"/*.frq; do
    sbatch << EOF
#!/bin/bash
#SBATCH --job-name=step5_fix_annot_${frq_file##*/}
#SBATCH --output=slurm_logs/step5_fix_annot_${frq_file##*/}_output.txt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=3
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=8:00:00
#SBATCH --get-user-env=L

# Run fix_annot script and output fixed frequency files to fixed_freq directory
python3 "$fix_annot_script_dir" -i "$frq_file" -o "$fixed_freq_dir/$(basename "${frq_file%.frq}")_fixed.frq"
EOF
done

# End of the script
