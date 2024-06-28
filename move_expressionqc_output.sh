#!/bin/bash

# Define the source folders
# "SRP064952" mist nog
folders=("SRP076426" "SRP077046" "SRP068609" "SRP063496" "ERP114636" "SRP129004" "SRP189239" "SRP155976" "SRP096757" "ERP109626" "SRP125961" "ERP113396")

# Define the base source path and destination path
base_source_path="/scratch/hb-functionalgenomics/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq"
destination_path="/scratch/hb-functionalgenomics/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/00-Final_files/expression_data"

# Copy the expressionqc_output folder from each source folder to the destination path
for folder in "${folders[@]}"; do
    source_path="$base_source_path/$folder/expressionqc_output"
    if [ -d "$source_path" ]; then
        dest_folder="$destination_path/$folder"
        mkdir -p "$dest_folder"
        cp -r "$source_path/." "$dest_folder/"
    else
        echo "$source_path does not exist."
    fi
done
