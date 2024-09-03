# Author: Peter Riesebos

import os
import subprocess

# Common parameters
VCF_FILE = "/scratch/hb-functionalgenomics/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/00-Final_files/final_sample_names_exp_filtered.vcf.gz"
EXPRESSION_DATA = "/scratch/hb-functionalgenomics/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/00-Final_files/merged_expression_data_zeros.txt.gz"
ANNOTATION_FILE = "/scratch/hb-functionalgenomics/projects/gut-bulk/reference/gencode_44_2023/annotation_file_build44_genes.tsv"
LINK_FILE = "/scratch/hb-functionalgenomics/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/00-Final_files/final_linkfile.txt"
MODE = "mbqtl"
PERMUTATIONS = 100
OUTPUT_PREFIX = "pub_rna_perm100_all"
MINGENOTYPECOUNT = 2

# Path to the jar file
JAR_PATH = "/scratch/hb-functionalgenomics/projects/gut-bulk/ongoing/users/umcg-priesebos/tools/MbQTL-1.5.0-SNAPSHOT-jar-with-dependencies.jar"

# Directory to store the generated sbatch scripts
OUTPUT_DIR = "mbqtl_sbatch_scripts"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# SBATCH script template
sbatch_script_template = """#!/bin/bash
#SBATCH --time=23:59:59
#SBATCH --mem=32g
#SBATCH --cpus-per-task=16
#SBATCH -J mbqtl_chr{chr}
#SBATCH -o mbqtl_logs/pub_rna_mbqtl_run_chr{chr}.log
#SBATCH -e mbqtl_logs/pub_rna_mbqtl_run_chr{chr}.err

module load java

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
        mingenotypecount=MINGENOTYPECOUNT
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
