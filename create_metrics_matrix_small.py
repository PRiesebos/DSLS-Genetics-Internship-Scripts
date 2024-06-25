import os
import gzip
import pandas as pd

# Define the base directory
base_dir = "/scratch/hb-functionalgenomics/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/5-expressieQC/all_alignment_output"

# Initialize a dictionary to store the data
data = {}

# Process multiple_metrics folder
multiple_metrics_dir = os.path.join(base_dir, "multiple_metrics")

for root, dirs, files in os.walk(multiple_metrics_dir):
    for file in files:
        if file == "multiple_metrics.alignment_summary_metrics.gz":
            sample = os.path.basename(root)
            file_path = os.path.join(root, file)
            
            with gzip.open(file_path, 'rt') as f:
                content = f.read()
                if "PAIR" in content:
                    for line in content.split("\n"):
                        if line.startswith(("PAIR", "UNPAIRED")):
                            values = line.split()
                            pct_pf_reads_aligned = values[6]  # PCT_PF_READS_ALIGNED is the 7th value in the line
                            if sample not in data:
                                data[sample] = {}
                            data[sample]['PCT_PF_READS_ALIGNED'] = pct_pf_reads_aligned

# Process rna_seq_metrics folder
rna_seq_metrics_dir = os.path.join(base_dir, "rna_seq_metrics")

for root, dirs, files in os.walk(rna_seq_metrics_dir):
    for file in files:
        if file.endswith("rnaseqmetrics.gz"):
            sample = os.path.basename(root)
            file_path = os.path.join(root, file)
            
            with gzip.open(file_path, 'rt') as f:
                lines = f.readlines()
                for i, line in enumerate(lines):
                    if line.strip() and not line.startswith("#"):
                        if i + 1 < len(lines) and lines[i + 1].strip() and not lines[i + 1].startswith("#"):
                            values = lines[i + 1].split()
                            pct_coding_bases = values[16]  # PCT_CODING_BASES is the 17th value in the row below the header
                            if sample not in data:
                                data[sample] = {}
                            data[sample]['PCT_CODING_BASES'] = pct_coding_bases
                            break

# Convert to a DataFrame for easier handling
df = pd.DataFrame.from_dict(data, orient='index').reset_index()
df.columns = ['Sample', 'PCT_PF_READS_ALIGNED', 'PCT_CODING_BASES']
print(df.head())

# Write the DataFrame to a gzipped text file
output_file = "/scratch/hb-functionalgenomics/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/5-expressieQC/all_alignment_output/alignment_metrics.txt.gz"
with gzip.open(output_file, 'wt') as f:
    df.to_csv(f, index=False, sep='\t')
