# Author: Peter Riesebos
# Purpose: Script used to submit mbQTL jobs

import os
import subprocess

# Common parameters
VCF_FILE = "/groups/umcg-fg/tmp04/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/datasets/combined/combined_vcf_files.vcf.gz"
EXPRESSION_DATA = "/groups/umcg-fg/tmp04/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/datasets/combined/combined_expression_matrix_protein_coding_filtered.txt.gz"
ANNOTATION_FILE = "/groups/umcg-fg/tmp04/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/references/gencode_44_2023/annotation_file_build44_genes_no_version.tsv"
LINK_FILE = "/groups/umcg-fg/tmp04/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/datasets/GTEx/tweaked_files/linkfile_gtex_final.txt"
MODE = "mbqtl"
PERMUTATIONS = 100
OUTPUT_PREFIX = "gtex_no_ver"
MINGENOTYPECOUNT = 2
MAF = 0.05

# Path to the jar file
JAR_PATH = "/groups/umcg-fg/tmp04/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/MbQTL-1.5.0-SNAPSHOT-jar-with-dependencies.jar"

# Directory to store the generated sbatch scripts
OUTPUT_DIR = "mbqtl_output_gtex_no_version"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# SBATCH script template
sbatch_script_template = """#!/bin/bash
#SBATCH --time=4:00:00
#SBATCH --mem=32gb
#SBATCH --cpus-per-task=16
#SBATCH -J mbqtl_logs/gtex_no_ver_{chr}
#SBATCH -o mbqtl_logs/gtex_no_ver_{chr}.log
#SBATCH -e mbqtl_logs/gtex_no_ver_{chr}.err

module load Java/11-LTS

java -jar {jar_path} \\
    -v {vcf_file} \\
    -e {expression_data} \\
    -a {annotation_file} \\
    -g {link_file} \\
    -m {mode} \\
    -perm {permutations} \\
    -o {output_prefix}_chr{chr} \\
    --mingenotypecount {mingenotypecount} \\
    --chr {chr} \\
    --maf {maf} \\
    --outputall

"""

# Generate and submit sbatch scripts for each chromosome
for chr_num in range(1, 23):  # Chromosomes 1 to 22
    sbatch_script = sbatch_script_template.format(
        chr=chr_num,
        jar_path=JAR_PATH,
        vcf_file=VCF_FILE,
        expression_data=EXPRESSION_DATA,
        annotation_file=ANNOTATION_FILE,
        link_file=LINK_FILE,
        mode=MODE,
        permutations=PERMUTATIONS,
        output_prefix=OUTPUT_PREFIX,
        mingenotypecount=MINGENOTYPECOUNT,
        maf=MAF
    )
    script_path = os.path.join(OUTPUT_DIR, f"submit_chr{chr_num}.sh")
    
    with open(script_path, "w") as f:
        f.write(sbatch_script)
    
    # Prepare the sbatch command
    sbatch_command = ["sbatch", script_path]
    
    # Submit the job using sbatch
    result = subprocess.run(sbatch_command, capture_output=True, text=True)
    output = result.stdout.strip()
    
    # Print the output for monitoring purposes
    print(output)

print("SLURM job submission scripts have been generated and submitted.")
