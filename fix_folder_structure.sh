#!/bin/bash

# Author: Peter Riesebos

# usage:
# ./create_folder_structure.sh /path/to/input/location
# input location example: /scratch/hb-functionalgenomics/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/SRP129004

# Check if the correct number of arguments were passed
if [ $# -ne 1 ]; then
  echo "Usage: $0 <input_location>"
  exit 1
fi

# Set the input location
input_location=$1

# Extract the last folder name from the input location path
root_folder_name=$(basename "$input_location")-test

# Derive the output location by removing the last folder name from the input location
output_location="${input_location%/*}"

# Create the root folder with the modified name in the output location
mkdir -p "$output_location/$root_folder_name"

# Get the list of samples (subfolders) in the rmats folder
samples=()
for d in "$input_location/rmats/"*; do
  if [ -d "$d" ]; then
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

# Copy FastQC files
for sample in "${samples[@]}"; do
  for file in "$input_location/$sample/$sample/"*; do
    cp "$file" "$output_location/$root_folder_name/fastqc/$sample/"
  done
done

# copy gvcf files
for file in "$input_location/gvcf/"*; do
  if [ -f "$file" ]; then
    sample_name=$(basename "$file" | cut -d '-' -f1)
    
    cp "$file" "$output_location/$root_folder_name/gvcf/$sample_name/"
  fi
done

# Copy leafcutter files
for sample in "${samples[@]}"; do
  for file in "$input_location/leafcutter/$sample/${sample}.junc.gz"; do
    cp "$file" "$output_location/$root_folder_name/leafcutter/$sample/"
  done
done

# Copy mark_duplicates files
for sample in "${samples[@]}"; do
  for file in "$input_location-logs/$sample/mark_duplicates/${sample}_duplicates.txt.gz"; do
    cp "$file" "$output_location/$root_folder_name/mark_duplicates/$sample/"
  done
done

# Copy multiple_metrics files
for sample in "${samples[@]}"; do
  for file in "$input_location-logs/$sample/$sample/"*; do
    if [ -f "$file" ]; then
      cp "$file" "$output_location/$root_folder_name/multiple_metrics/$sample/"
    fi
  done
done

# Copy rmats files
for sample in "${samples[@]}"; do
  for file in "$input_location/rmats/$sample/"*; do
    if [ -f "$file" ]; then
      cp "$file" "$output_location/$root_folder_name/rmats/$sample/"
    fi
  done
done

# copy rna_seq_metrics files
for file in "$input_location/rna_seq_metrics/"*; do
  if [ -f "$file" ]; then
    sample_name=$(basename "$file" | sed 's/\..*//' | cut -d '_' -f1)

    cp "$file" "$output_location/$root_folder_name/rna_seq_metrics/$sample_name/"
  fi
done

# Copy star files
for sample in "${samples[@]}"; do
  for file in "$input_location-logs/$sample/star/"*; do
    if [ -f "$file" ]; then
      cp "$file" "$output_location/$root_folder_name/star/$sample/"
    fi
  done
done

echo "All done! Hasta La Vista."
