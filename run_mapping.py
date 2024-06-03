# Author: Peter Riesebos

import os
import subprocess

# Define constants
sbatch_script_template = """#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=6:00:00
#SBATCH --mem=2g
#SBATCH --cpus-per-task=2
#SBATCH -J JOB_chr{chromosome}
#SBATCH -o sbatch_scripts/chr{chromosome}.log
#SBATCH -e sbatch_scripts/chr{chromosome}.err

set -e
set -u

ml Java
threads=2
nice -n20 java -Xmx1500m \\
        -Djava.util.concurrent.ForkJoinPool.common.parallelism=$threads \\
        -Dmaximum.threads=$threads -Dthread.pool.size=$threads \\
        -jar /scratch/hb-functionalgenomics/projects/gut-bulk/ongoing/users/umcg-priesebos/tools/MbQTL-1.4.2-SNAPSHOT-jar-with-dependencies.jar \\
        -v /scratch/hb-functionalgenomics/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/6-imputatie/output/postimpute/chr{chromosome}.dose.filtered.vcf.gz \\
        -e /scratch/hb-functionalgenomics/projects/gut-bulk/reference/expression_data/combinedHealthyTissue_TPM_Log2_QQ_CovCor_Exp_No_Version.txt.gz \\
        -a /scratch/hb-functionalgenomics/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/6-imputatie/output/postimpute/annotation_file.tsv \\
        -g /scratch/hb-functionalgenomics/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/6-imputatie/output/postimpute/linkfile.txt \\
        -m mbqtl \\
        --perm 10 \\
        --chr {chromosome} \\
        -o chr{chromosome}.mapped
"""

# Define the range of chromosomes
chromosome_range = range(1, 23)

# Directory to store the generated sbatch scripts
output_dir = "sbatch_scripts"
os.makedirs(output_dir, exist_ok=True)

# Variable to store the job ID of the previous job
previous_job_id = None

# Generate and submit sbatch scripts
for chromosome in chromosome_range:
    sbatch_script = sbatch_script_template.format(chromosome=chromosome)
    script_path = os.path.join(output_dir, f"submit_chr{chromosome}.sh")
    
    with open(script_path, "w") as f:
        f.write(sbatch_script)
    
    # Prepare the sbatch command
    sbatch_command = ["sbatch"]
    if previous_job_id is not None:
        sbatch_command.append(f"--dependency=afterok:{previous_job_id}")
    sbatch_command.append(script_path)
    
    # Submit the job using sbatch and capture the job ID
    result = subprocess.run(sbatch_command, capture_output=True, text=True)
    output = result.stdout.strip()
    
    # Extract the job ID from the output
    job_id = output.split()[-1]
    previous_job_id = job_id

print("SLURM job submission scripts have been generated and submitted.")
