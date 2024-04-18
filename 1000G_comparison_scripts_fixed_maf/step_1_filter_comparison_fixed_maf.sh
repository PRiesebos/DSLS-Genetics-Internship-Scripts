#!/bin/bash

# Define variables
CRS=(0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0)
DPS=(1 2 3 4 5 6 7 8 9 10)

output_dir='/scratch/hb-functionalgenomics/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/SRP129004/filter_comparison_fixed_maf'
output_vcf_dir="$output_dir/vcf"
vcf_file_dir='/scratch/hb-functionalgenomics/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/SRP129004/genotyping_output/split_merged_output_dotrevive_annot.vcf.gz'
filter_script_dir='/scratch/hb-functionalgenomics/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/extra_scripts/custom_vcf_filter.py'

# Step 1: Run filter script with varying CR parameters
for cr in "${CRS[@]}"; do
    # Step 2: Run filter script with varying DP parameters
    for dp in "${DPS[@]}"; do
        sbatch << EOF
#!/bin/bash
#SBATCH --job-name=step1_cr_${cr}_dp_${dp}
#SBATCH --output=slurm_logs/step1_cr_${cr}_dp_${dp}_output.txt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=23:59:59
#SBATCH --get-user-env=L

python3 "$filter_script_dir" -i "$vcf_file_dir" -o "$output_vcf_dir/cr_${cr}_dp_${dp}_output.vcf.gz" -cr "$cr" -maf 0.0 -dp "$dp" -gq 10 --no_snv_vqsr_check --no_indel_vqsr_check --remove_non_pass_snv --remove_non_pass_indel --replace_poor_quality_genotypes
EOF
    done
done

# End of the script
