import os
import argparse
import re

def main(input_folder):
    # Define file paths based on the input folder
    combined_file = os.path.join(input_folder, 'combined_samples_PCs.txt')
    mapping_file = os.path.join(input_folder, 'study_samples_mapping.txt')
    output_file = os.path.join(input_folder, 'combined_samples_PCs_with_study.txt')

    # Step 1: Read combined_samples_PCs.txt
    combined_data = []

    with open(combined_file, 'r') as f:
        header = f.readline().strip().split('\t')
        for line in f:
            parts = line.strip().split('\t')
            combined_data.append(parts)

    # Step 2: Read study_samples_mapping.txt
    sample_to_study = {}

    with open(mapping_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            study_name = parts[0]
            samples = parts[1:]
            for sample in samples:
                # Extract the part of the sample name that matches SRR or ERR followed by numbers
                match = re.match(r'(SRR|ERR)\d+', sample)
                if match:
                    sample_name = match.group(0)
                    sample_to_study[sample_name] = study_name

    # Step 3: Modify combined_samples_PCs.txt data with study names
    modified_data = []
    modified_header = header.copy()
    modified_header.insert(0, 'Study')

    for row in combined_data:
        sample_name = row[0]
        match = re.match(r'(SRR|ERR)\d+', sample_name)
        if match:
            sample_key = match.group(0)
            if sample_key in sample_to_study:
                study_name = sample_to_study[sample_key]
                modified_row = [study_name] + row
                modified_data.append(modified_row)
            else:
                # Handle cases where sample name prefix does not have a corresponding study name
                print(f"Sample '{sample_name}' not found in study mapping.")
        else:
            # Handle cases where sample name does not match expected pattern
            print(f"Sample '{sample_name}' does not match expected pattern (SRR or ERR followed by numbers).")

    # Step 4: Write modified data to a new file
    with open(output_file, 'w') as f:
        # Write header
        f.write('\t'.join(modified_header) + '\n')
        # Write data rows
        for row in modified_data:
            f.write('\t'.join(map(str, row)) + '\n')

    print(f"Modified data has been written to '{output_file}'.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Add study names to combined samples PCs file based on mapping.')
    parser.add_argument('input_folder', type=str, help='Folder path containing input files')
    args = parser.parse_args()

    main(args.input_folder)
