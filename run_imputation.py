import os
import subprocess

# Define the SLURM script template
sbatch_script_template = """#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=6G
#SBATCH --job-name="ImputeGenotypes_{study_id}"
#SBATCH -o sbatch_scripts/imputation/{study_id}.log
#SBATCH -e sbatch_scripts/imputation/{study_id}.err

module load Java/17.0.6

export SINGULARITY_CACHEDIR=/scratch/hb-functionalgenomics/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/6-imputatie/eQTLGenImpute/cache/singularity
export NXF_HOME=/scratch/hb-functionalgenomics/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/6-imputatie/eQTLGenImpute/cache/nextflow

nextflow_path=/scratch/hb-functionalgenomics/projects/gut-bulk/ongoing/users/umcg-priesebos/tools
reference_path=/scratch/hb-functionalgenomics/projects/gut-bulk/reference/eQTLGenP2/hg38

cohort_name={cohort_name}
qc_input_folder={qc_input_folder}
output_path={output_path}
genome_build="GRCh38"

NXF_VER=21.10.6 ${{nextflow_path}}/nextflow run /scratch/hb-functionalgenomics/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/6-imputatie/eQTLGenImpute/eQTLGenImpute.nf \\
--qcdata ${{qc_input_folder}} \\
--target_ref ${{reference_path}}/genome_reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa \\
--ref_panel_hg38 ${{reference_path}}/harmonizing_reference/30x-GRCh38_NoSamplesSorted \\
--eagle_genetic_map ${{reference_path}}/phasing_reference/genetic_map/genetic_map_hg38_withX.txt.gz \\
--eagle_phasing_reference ${{reference_path}}/phasing_reference/phasing/ \\
--minimac_imputation_reference ${{reference_path}}/imputation_reference/ \\
--cohort_name ${{cohort_name}} \\
--genome_build ${{genome_build}} \\
--outdir ${{output_path}} \\
-profile slurm,singularity \\
-resume
"""

# Define the list of study IDs
# moet nog     "SRP113470/genotypes/final_beagle/",
study_ids = [
    "ERP114636/genotypes/final_beagle/",
    "SRP189239/genotypes/final_beagle/",
]

base_qc_input_folder = "/scratch/hb-functionalgenomics/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq"
base_output_path = "/scratch/hb-functionalgenomics/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq"

# Directory to store the generated sbatch scripts
output_dir = "sbatch_scripts/imputation"
os.makedirs(output_dir, exist_ok=True)

# Variable to store the job ID of the previous job
previous_job_id = None

# Generate and submit sbatch scripts
for study_id in study_ids:
    main_study_id = study_id.split('/')[0]  # Get the first part of the study_id
    cohort_name = main_study_id  # Use the main_study_id as the cohort name
    qc_input_folder = f"{base_qc_input_folder}/{study_id}"
    output_path = f"{base_output_path}/{main_study_id}/imputation_output/"

    # Create output directory if it doesn't exist
    os.makedirs(output_path, exist_ok=True)

    # Generate the sbatch script content
    sbatch_script = sbatch_script_template.format(
        study_id=main_study_id,
        cohort_name=cohort_name,
        qc_input_folder=qc_input_folder,
        output_path=output_path
    )
    script_path = os.path.join(output_dir, f"submit_{main_study_id}.sh")
    
    # Write the sbatch script to a file
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
    
    print(f"Submitted job for {main_study_id} with job ID {job_id}")

print("SLURM job submission scripts have been generated and submitted.")