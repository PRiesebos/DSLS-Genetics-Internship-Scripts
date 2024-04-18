#!/bin/bash

# Define variables
MAFS=(0.01 0.007)
CRS=(0.1 0.2)
DPS=(1 2)
output_dir='/scratch/hb-functionalgenomics/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/SRP129004/filter_comparison_test'
output_vcf_dir="$output_dir/vcf"
vcf_file_dir='/scratch/hb-functionalgenomics/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/SRP129004/genotyping_output/split_merged_output_dotrevive_annot.vcf.gz'
filter_script_dir='/scratch/hb-functionalgenomics/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/extra_scripts/custom_vcf_filter.py'

# Step 1: Run filter script with varying MAF parameters
for maf in "${MAFS[@]}"; do
    sbatch << EOF
#!/bin/bash
#SBATCH --job-name=step1_maf_${maf}
#SBATCH --output=slurm_logs/step1_maf_${maf}_output.txt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=8:00:00
#SBATCH --get-user-env=L

python3 "$filter_script_dir" -i "$vcf_file_dir" -o "$output_vcf_dir/maf_${maf}_output.vcf.gz" -cr 0.5 -maf "$maf" -dp 5 -gq 10 --no_snv_vqsr_check --no_indel_vqsr_check --remove_non_pass_snv --remove_non_pass_indel --replace_poor_quality_genotypes
EOF
done &

# Step 2: Run filter script with varying CR parameters
for cr in "${CRS[@]}"; do
    sbatch << EOF
#!/bin/bash
#SBATCH --job-name=step2_cr_${cr}
#SBATCH --output=slurm_logs/step2_cr_${cr}_output.txt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=8:00:00
#SBATCH --get-user-env=L

python3 "$filter_script_dir" -i "$vcf_file_dir" -o "$output_vcf_dir/cr_${cr}_output.vcf.gz" -cr "$cr" -maf 0.01 -dp 5 -gq 10 --no_snv_vqsr_check --no_indel_vqsr_check --remove_non_pass_snv --remove_non_pass_indel --replace_poor_quality_genotypes
EOF
done &

# Step 3: Run filter script with varying DP parameters
for dp in "${DPS[@]}"; do
    sbatch << EOF
#!/bin/bash
#SBATCH --job-name=step3_dp_${dp}
#SBATCH --output=slurm_logs/step3_dp_${dp}_output.txt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=8:00:00
#SBATCH --get-user-env=L

python3 "$filter_script_dir" -i "$vcf_file_dir" -o "$output_vcf_dir/dp_${dp}_output.vcf.gz" -cr 0.5 -maf 0.01 -dp "$dp" -gq 10 --no_snv_vqsr_check --no_indel_vqsr_check --remove_non_pass_snv --remove_non_pass_indel --replace_poor_quality_genotypes
EOF
done

# End of the script
