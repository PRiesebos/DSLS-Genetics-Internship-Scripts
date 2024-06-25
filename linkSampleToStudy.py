import os
import glob

# Directory where the studies are located
root_dir = '/scratch/hb-functionalgenomics/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/'

# Output file to store study names and sample names
output_file = 'study_samples_mapping.txt'

# Open the output file in write mode
with open(output_file, 'w') as f:
    # Iterate over each study directory
    for study_dir in os.listdir(root_dir):
        study_path = os.path.join(root_dir, study_dir)
        
        # Check if it's a directory and start with 'SRP' or 'ERP'
        if os.path.isdir(study_path) and (study_dir.startswith('SRP') or study_dir.startswith('ERP')):
            # Initialize lists to store sample names
            srr_samples = []
            err_samples = []
            
            # Search for *ReadsPerGene.out.tab.gz files recursively within the study directory
            for root, dirs, files in os.walk(study_path):
                for file in files:
                    if file.endswith('ReadsPerGene.out.tab.gz'):
                        file_path = os.path.join(root, file)
                        
                        # Determine if the sample starts with 'SRR' or 'ERR'
                        sample_name = os.path.basename(file).split('ReadsPerGene.out.tab.gz')[0]
                        if sample_name.startswith('SRR'):
                            srr_samples.append(sample_name)
                        elif sample_name.startswith('ERR'):
                            err_samples.append(sample_name)
            
            # Write study name and sample names to output file
            if srr_samples:
                f.write(f"{study_dir}\t{' '.join(srr_samples)}\n")
            if err_samples:
                f.write(f"{study_dir}\t{' '.join(err_samples)}\n")

print(f"Study names and sample names written to '{output_file}'.")
