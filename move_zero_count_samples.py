import os
import shutil

# Directory paths
star_dir = "/scratch/hb-functionalgenomics/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/5-expressieQC/all_alignment_output/star"
multiple_metrics_dir = "/scratch/hb-functionalgenomics/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/5-expressieQC/all_alignment_output/multiple_metrics"
rna_seq_metrics_dir = "/scratch/hb-functionalgenomics/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/5-expressieQC/all_alignment_output/rna_seq_metrics"

# Directory to move removed samples
removed_samples_parent_dir = "/scratch/hb-functionalgenomics/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/5-expressieQC/removed_samples"

# Samples to be removed
samples_to_remove = [
    "ERR3262444", "ERR3262424", "ERR3079564", "ERR3262401", "ERR3262421", "ERR3079595",
    "ERR3262416", "ERR3079567", "ERR3079584", "ERR3262431", "ERR3262423", "ERR3079596",
    "ERR3079570", "ERR3262419", "ERR3262435", "ERR3079586", "ERR3262429", "ERR3079592",
    "ERR3262415", "ERR3079572", "ERR3262430", "ERR3262405", "ERR3079578", "ERR3262408",
    "ERR3262432", "ERR3262404", "ERR3262441", "ERR3079587", "ERR3079582", "ERR3079577",
    "ERR3079568", "ERR3079561", "ERR3262427", "ERR3079560", "ERR3262411", "ERR3262402",
    "ERR3262407", "ERR3079579", "ERR3079569", "ERR3262440", "ERR3262413", "ERR3262414",
    "ERR3079580", "ERR3262434", "ERR3262442", "ERR3079565", "ERR3079585", "ERR3262433",
    "ERR3262436", "ERR3079599", "ERR3079594", "ERR3079597", "ERR3079563", "ERR3079559",
    "ERR3262417", "ERR3262426", "ERR3262400", "ERR3079573", "ERR3079588", "ERR3262412",
    "ERR3262420", "ERR3079575", "ERR3262399", "ERR3079589", "ERR3262425", "ERR3262443",
    "ERR3262428", "ERR3079576", "ERR3079583"
]

# Function to check if a folder name corresponds to a sample that should be removed
def should_remove(folder_name):
    for sample in samples_to_remove:
        if sample in folder_name:
            return True
    return False

# Function to move removed samples to removed_samples_parent_dir
def move_removed_samples(root_dir, sub_dir):
    print(f"Moving removed samples from {root_dir} to {removed_samples_parent_dir}/{sub_dir}")
    for folder_name in os.listdir(root_dir):
        full_path = os.path.join(root_dir, folder_name)
        if os.path.isdir(full_path) and should_remove(folder_name):
            print(f"Moving {folder_name}")
            destination_path = os.path.join(removed_samples_parent_dir, sub_dir, folder_name)
            os.makedirs(destination_path, exist_ok=True)  # Create parent directories if needed
            shutil.move(full_path, destination_path)

# Move removed samples from each directory to their respective subdirectories in removed_samples_parent_dir
move_removed_samples(star_dir, "star")
move_removed_samples(multiple_metrics_dir, "multiple_metrics")
move_removed_samples(rna_seq_metrics_dir, "rna_seq_metrics")

print("Done")
