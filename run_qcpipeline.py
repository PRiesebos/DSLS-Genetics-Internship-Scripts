# Author: Peter Riesebos

import os
import subprocess

# List of input directories
input_dirs = ["SRP068609","SRP064952","SRP063496", "SRP076426", "SRP077046", "SRP096757", "SRP113470", "SRP189239", "SRP048801"]

# Nextflow pipeline script
nf_script = "/scratch/hb-functionalgenomics/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/QCPipeline/main.nf"

# SBATCH script template
sbatch_script_template = """#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=12:00:00
#SBATCH --mem=2g
#SBATCH --cpus-per-task=1
#SBATCH -J QC_{input_dir}
#SBATCH -o sbatch_scripts/QC/{input_dir}.log
#SBATCH -e sbatch_scripts/QC/{input_dir}.err

set -e
set -u

module load Java

# Create the output directory
output_dir="/scratch/hb-functionalgenomics/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/{input_dir}-build44/genotypes"
mkdir -p "$output_dir"

# Create a temporary nextflow.config file with the current input directory
temp_config="temp_nextflow_{input_dir}.config"
cat <<EOT > $temp_config
nextflow.enable.dsl=2

process.executor = "slurm"
process.container = "/scratch/hb-functionalgenomics/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/QCPipeline/pub-rna-nocache_latest.sif"

// Adjust the params to your usecase.
params {{
    inputDir = "/scratch/hb-functionalgenomics/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/{input_dir}-build44/genotypes"
    maf = 0.01 
    geno = 0.05 
    mind = 0.05 
    hwe = 1e-6
    kingTableFilter = 0.04419417382
    populationOutlierThreshold = 1
    refPath = '/scratch/hb-functionalgenomics/projects/gut-bulk/reference/1000g_HJ/concatenated_all.bed'
    refPop = '/scratch/hb-functionalgenomics/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/QCPipeline/bin/ref_pop.txt'
    refAFs = '/scratch/hb-functionalgenomics/projects/gut-bulk/reference/1000g_HJ/concatenated_all.afreq'
    annotationGTF = '/scratch/hb-functionalgenomics/projects/gut-bulk/reference/gencode_44_2023/gencode.v44.primary_assembly.annotation.gtf'
    beagleJarDir = '/scratch/hb-functionalgenomics/projects/gut-bulk/ongoing/users/umcg-priesebos/tools/beagle.27Jan18.7e1.jar'
    mapFile = '/scratch/hb-functionalgenomics/projects/gut-bulk/reference/plink_map/combined.map'
    callRate = 0.5
}}

singularity {{
    enabled = true
    autoMounts = true
    runOptions = "--bind /scratch/"
    cacheDir = "/scratch/hb-functionalgenomics/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/combined_pipeline/apptainer_cache"
}}
EOT

# Run the Nextflow job
nextflow run -c $temp_config {nf_script}
"""

# Directory to store the generated sbatch scripts
output_dir = "sbatch_scripts"
os.makedirs(output_dir, exist_ok=True)

# Variable to store the job ID of the previous job
previous_job_id = None

# Generate and submit sbatch scripts
for input_dir in input_dirs:
    sbatch_script = sbatch_script_template.format(input_dir=input_dir, nf_script=nf_script)
    script_path = os.path.join(output_dir, f"submit_{input_dir}.sh")
    
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
