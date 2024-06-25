#!/bin/bash

# Author: Peter Riesebos

# usage:
# ./create_folder_structure.sh /path/to/input/location
# input location example: /scratch/hb-functionalgenomics/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/SRP129004
# no / at the end of the input path!!

# Check if the correct number of arguments were passed
if [ $# -ne 1 ]; then
  echo "Usage: $0 <input_location>"
  exit 1
fi

# Set the input location and remove any trailing slash
input_location=${1%/}

# Extract the last folder name from the input location path
root_folder_name=$(basename "$input_location")

# Remove the '-build44' suffix from the folder name if it exists
root_folder_name=${root_folder_name%-build44}

# Derive the output location by removing the last folder name from the input location
output_location="${input_location%/*}"

# Create the root folder with the modified name in the output location
mkdir -p "$output_location/$root_folder_name"

# Get the list of samples (subfolders) in the rmats folder
samples=()
for d in "$input_location/"*; do
  if [ -d "$d" ] && [ "${d##*/}" != "genotypes" ]; then
    samples+=("${d##*/}")
  fi
done

# List of subfolder names
subfolder_names=("fastqc" "gvcf" "leafcutter" "mark_duplicates" "multiple_metrics" "rmats" "rna_seq_metrics" "star")

# Loop through the list of subfolder names
for subfolder_name in "${subfolder_names[@]}"; do
  # Create the subfolder inside the root folder
  mkdir -p "$output_location/$root_folder_name/$subfolder_name"
  
  # Loop through the samples
  for sample in "${samples[@]}"; do
    # Create the sample folder inside the subfolder
    mkdir -p "$output_location/$root_folder_name/$subfolder_name/$sample"
    echo "Created folder: $output_location/$root_folder_name/$subfolder_name/$sample"
  done
done

echo "1/10. copying FastQC files"
# Copy FastQC files
for sample in "${samples[@]}"; do
  for file in "$input_location/$sample/fastqc/"*; do
    cp "$file" "$output_location/$root_folder_name/fastqc/$sample/"
  done
done

echo "2/10. copying gvcf files"
# copy gvcf files
for sample in "${samples[@]}"; do
  for file in "$input_location/$sample/gvcf/"*; do
    cp "$file" "$output_location/$root_folder_name/gvcf/$sample/"
  done
done

echo "3/10. moving leafcutter files"
# Move leafcutter files
for sample in "${samples[@]}"; do
  for file in "$input_location/$sample/$sample/${sample}.junc.gz"; do
    mv "$file" "$output_location/$root_folder_name/leafcutter/$sample/"
  done
done

echo "4/10. copying mark_duplicates files"
# Copy mark_duplicates files
for sample in "${samples[@]}"; do
  for file in "$input_location/$sample/mark_duplicates/${sample}_duplicates.txt.gz"; do
    cp "$file" "$output_location/$root_folder_name/mark_duplicates/$sample/"
  done
done

echo "5/10. moving multiple_metrics files"
# Move multiple_metrics files
for sample in "${samples[@]}"; do
  for file in "$input_location/$sample/$sample/multiple_metrics."*; do
    if [ -f "$file" ]; then
      mv "$file" "$output_location/$root_folder_name/multiple_metrics/$sample/"
    fi
  done
done

echo "6/10. copying rmats files"
# Copy rmats files
for sample in "${samples[@]}"; do
  for file in "$input_location/$sample/$sample/"*; do
    if [ -f "$file" ]; then
      cp "$file" "$output_location/$root_folder_name/rmats/$sample/"
    fi
  done
done

echo "7/10. copying rna_seq_metrics files"
# Copy rna_seq_metrics files
for sample in "${samples[@]}"; do
  for file in "$input_location/$sample/${sample}_rnaseqmetrics.gz"; do
    cp "$file" "$output_location/$root_folder_name/rna_seq_metrics/$sample/"
  done
done

echo "8/10. copying star files"
# Copy star files
for sample in "${samples[@]}"; do
  for file in "$input_location/$sample/star/"*; do
    cp "$file" "$output_location/$root_folder_name/star/$sample/"
  done
done

echo "9/10. copying genotypes files"
# Copy genotypes folder
mkdir -p "$output_location/$root_folder_name/genotypes"
cp -r "$input_location/genotypes/"* "$output_location/$root_folder_name/genotypes/"

echo "10/10. All done! Hasta La Vista."
