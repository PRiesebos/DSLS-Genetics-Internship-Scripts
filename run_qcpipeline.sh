#!/bin/bash

# List of input directories
input_dirs=("SRP64952" "SRP076426" "SRP077046" "SRP096757" "SRP113470" "genotyping_output/SRP125961" "genotyping_output/SRP155976" "SRP189239")

# Template nextflow.config file
template_config="/scratch/hb-functionalgenomics/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/QCPipeline/nextflow.config"

# Nextflow pipeline script
nf_script="/scratch/hb-functionalgenomics/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/QCPipeline/main.nf"

for input_dir in "${input_dirs[@]}"; do
  # Create a temporary nextflow.config file with the current input directory
  cat $template_config | sed "s|inputDir = .*|inputDir = \"/scratch/hb-functionalgenomics/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/${input_dir}\"|g" > temp_nextflow.config

  # Run the Nextflow job
  nextflow run -c temp_nextflow.config $nf_script

  # Sleep for 1 hour
  sleep 1h
done