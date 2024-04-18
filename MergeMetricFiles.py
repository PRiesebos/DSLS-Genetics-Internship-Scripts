"""
Author: Peter Riesebos
"""

import os
import argparse
import gzip
import csv

def merge_files(directory_path):
    # Get a list of all files in the directory ending with '_outfile.txt.gz'
    files = [file for file in os.listdir(directory_path) if file.endswith('_outfile.txt.gz')]

    # Initialize an empty dictionary to store the data
    data_dict = {}

    # Iterate through each file
    for file in files:
        file_path = os.path.join(directory_path, file)
        with gzip.open(file_path, 'rt', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                metric = row.pop('-')  # Remove the '-' column and store its value in 'metric'
                if metric not in data_dict:
                    data_dict[metric] = {}
                data_dict[metric].update(row)

    # Create a list of dictionaries for each row in the final result
    rows = [{'metric': metric, **data} for metric, data in data_dict.items()]

    # Define the output file path
    output_path = os.path.join(directory_path, 'metric_files_merged.tsv.gz')

    # Save the merged data to a gzipped TSV file
    with gzip.open(output_path, 'wt', encoding='utf-8', newline='') as f:
        fieldnames = list(data_dict[next(iter(data_dict))].keys())
        writer = csv.DictWriter(f, fieldnames=['metric'] + fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)

    print(f"Combined data saved to {output_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Merge files in a specified directory.')
    parser.add_argument('-i', '--input_path', type=str, help='Path to the folder with files')

    args = parser.parse_args()

    if args.input_path:
        merge_files(args.input_path)
    else:
        print("Please provide the input path using the -i or --input_path argument.")
