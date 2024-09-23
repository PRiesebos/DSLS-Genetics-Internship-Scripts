#!/bin/bash
#SBATCH --job-name=vcf_merge
#SBATCH --output=vcf_merge_%j.out
#SBATCH --error=vcf_merge_%j.err
#SBATCH --time=23:59:59   # Set an appropriate time limit for your job
#SBATCH --mem=48gb          # Set the memory requirement for your job
#SBATCH --cpus-per-task=16 # Number of CPU cores to allocate
#SBATCH --nodes=1

# Author: Peter Riesebos
# Purpose: Script used for submitting an sbatch job to SLURM for merging several vcf files.

# Load the BCFtools module (adjust based on your environment)
module load BCFtools

# Define input and output files
VCF1="/groups/umcg-fg/tmp04/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/datasets/pub_rna/final_files_pub_rna/final_sample_names_exp_filtered_30_cutoff.vcf.gz"
VCF2="/groups/umcg-fg/tmp04/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/datasets/GTEx/2024-06-24-WGS_Imputed/sorted_merged_output.vcf.gz"
VCF3="/groups/umcg-fg/tmp04/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/datasets/Werna/genotype/werna_merged_filtered_chrs.vcf.gz"
OUTPUT="/groups/umcg-fg/tmp04/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/datasets/combined/combined_vcf_files.vcf.gz"

# Merge the VCF files using BCFtools
bcftools merge $VCF1 $VCF2 $VCF3 -Oz -o $OUTPUT --threads 16

# Print a message to indicate completion
echo "VCF files merged successfully into $OUTPUT"