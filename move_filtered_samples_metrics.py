import os
import shutil
import sys

# Directory paths
def get_directory_paths(base_dir):
    return {
        "star": os.path.join(base_dir, "star"),
        "multiple_metrics": os.path.join(base_dir, "multiple_metrics"),
        "rna_seq_metrics": os.path.join(base_dir, "rna_seq_metrics"),
        "removed_samples_parent_dir": os.path.join(base_dir, "removed_samples")
    }

# Samples to be removed
samples_to_remove = [
    "SRR8774204", "SRR8774209", "SRR8774211"
]

# Function to check if a folder name corresponds to a sample that should be removed
def should_remove(folder_name):
    return any(sample in folder_name for sample in samples_to_remove)

# Function to move removed samples to removed_samples_parent_dir
def move_removed_samples(root_dir, sub_dir, removed_samples_parent_dir):
    print(f"Moving removed samples from {root_dir} to {removed_samples_parent_dir}/{sub_dir}")
    for folder_name in os.listdir(root_dir):
        full_path = os.path.join(root_dir, folder_name)
        if os.path.isdir(full_path) and should_remove(folder_name):
            print(f"Moving {folder_name}")
            destination_path = os.path.join(removed_samples_parent_dir, sub_dir, folder_name)
            os.makedirs(destination_path, exist_ok=True)  # Create parent directories if needed
            shutil.move(full_path, destination_path)

# Function to create necessary directories if they do not exist
def create_directories(directories):
    for key, path in directories.items():
        if key != "removed_samples_parent_dir":  # Skip the main parent dir for individual creation
            os.makedirs(path, exist_ok=True)
    # Create subdirectories for removed samples
    os.makedirs(os.path.join(directories["removed_samples_parent_dir"], "star"), exist_ok=True)
    os.makedirs(os.path.join(directories["removed_samples_parent_dir"], "multiple_metrics"), exist_ok=True)
    os.makedirs(os.path.join(directories["removed_samples_parent_dir"], "rna_seq_metrics"), exist_ok=True)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <base_dir>")
        sys.exit(1)

    base_dir = sys.argv[1]
    directories = get_directory_paths(base_dir)
    
    # Create directories if they do not exist
    create_directories(directories)
    
    # Move removed samples from each directory to their respective subdirectories in removed_samples_parent_dir
    move_removed_samples(directories["star"], "star", directories["removed_samples_parent_dir"])
    move_removed_samples(directories["multiple_metrics"], "multiple_metrics", directories["removed_samples_parent_dir"])
    move_removed_samples(directories["rna_seq_metrics"], "rna_seq_metrics", directories["removed_samples_parent_dir"])

    print("Done")
