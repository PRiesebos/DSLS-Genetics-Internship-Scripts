import os

# Directory paths
star_dir = "/scratch/hb-functionalgenomics/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/5-expressieQC/all_alignment_output/star"
multiple_metrics_dir = "/scratch/hb-functionalgenomics/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/5-expressieQC/all_alignment_output/multiple_metrics"
rna_seq_metrics_dir = "/scratch/hb-functionalgenomics/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/5-expressieQC/all_alignment_output/rna_seq_metrics"


# Read first_parts.txt into a set
keep_samples = set()
with open("/scratch/hb-functionalgenomics/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/5-expressieQC/first_parts.txt", "r") as f:
    for line in f:
        keep_samples.add(line.strip())

# Function to check if a folder name should be kept
def should_keep(folder_name):
    for sample in keep_samples:
        if sample in folder_name:
            return True
    return False

# Function to remove unwanted folders
def remove_unwanted_folders(root_dir):
    print(f"Removing unwanted folders in {root_dir}")
    for folder_name in os.listdir(root_dir):
        if os.path.isdir(os.path.join(root_dir, folder_name)) and not should_keep(folder_name):
            print(f"Removing {folder_name}")
            os.system(f"rm -rf '{os.path.join(root_dir, folder_name)}'")

# Remove unwanted folders in each directory
remove_unwanted_folders(star_dir)
remove_unwanted_folders(multiple_metrics_dir)
remove_unwanted_folders(rna_seq_metrics_dir)

print("Done")
