#!/bin/bash

# List of identifiers
identifiers=("SRP063496" "SRP064952" "SRP068609" "SRP076426" "SRP077046" "SRP096757" "SRP113470" "SRP189239")

# Loop over each identifier
for id in "${identifiers[@]}"; do
  # Construct the directory and tar names
  dir_name="${id}-test"
  tar_name="${id}-logs.tar.gz"
  
  # Run the tar command
  tar -zvcf "$tar_name" "$dir_name"
done
