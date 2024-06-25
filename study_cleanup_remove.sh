#!/bin/bash

# List of identifiers
# deze moeten nog:  "SRP077046" "SRP096757" "SRP113470" "SRP189239"
identifiers=("SRP077046" "SRP096757" "SRP113470" "SRP189239")

# Loop over each identifier
for id in "${identifiers[@]}"; do
  # Construct the directory and tar names
  dir_name="${id}-test"
  
  # Run the tar command
  rm -rf "$dir_name"
done
