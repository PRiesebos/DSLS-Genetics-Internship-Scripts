#!/bin/bash

# List of input directories
input_dirs=("SRP063496" "SRP068609" "SRP076426" "SRP077046" "SRP096757" "SRP113470" "SRP189239")

# Template nextflow.config file
template_config="/scratch/hb-functionalgenomics/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/combined_pipeline/pub-rna/nextflow.config"

# Nextflow pipeline script
nf_script="/scratch/hb-functionalgenomics/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/combined_pipeline/pub-rna/main.nf"

for input_dir in "${input_dirs[@]}"; do
  # Create the output directory
  output_dir="/scratch/hb-functionalgenomics/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/${input_dir}-build44"
  mkdir -p "$output_dir"

  # Create a temporary nextflow.config file with the current input directory
  temp_config="temp_nextflow_${input_dir}.config"
  sed -e "s|out_dir = .*|out_dir = \"$output_dir\"|g" \
      -e "s|input_txt = .*|input_txt = \"/scratch/hb-functionalgenomics/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/combined_pipeline/pub-rna/test_data/rna_pub_${input_dir}.txt\"|g" \
      $template_config > $temp_config

  # Run the Nextflow job in the background
  nextflow run -c $temp_config $nf_script &
  
  sleep 2h

  # Optionally remove the temporary config file if no longer needed
  # rm -f $temp_config
done
