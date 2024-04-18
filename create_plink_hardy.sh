#!/bin/bash

# Create output directory if it doesn't exist
output_dir="../hardy"
mkdir -p "$output_dir"

# Iterate over each .vcf.gz file in the current directory
for vcf_file in *filtered.vcf.gz; do
    # Extract filename without extension
    filename=$(basename -- "$vcf_file")
    filename_no_ext="${filename%.*}"

    # Run PLINK 2 missingness analysis
    plink2 --vcf "$vcf_file" --hardy --out "$output_dir/$filename_no_ext" --allow-misleading-out-arg
done
