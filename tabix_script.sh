#!/bin/bash

# Set the common part of the filenames
common_part="chr"

# Loop through all files with the common part and run tabix
for file in "${common_part}"*; do
    # Check if the file exists
    if [ -e "${file}" ]; then
        # Run tabix on the file
        tabix -p vcf "${file}"
    else
        echo "File ${file} not found."
    fi
done
