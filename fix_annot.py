"""
Author: Peter Riesebos
"""

import argparse

def process_line(line):
    parts = line.strip().split()
    if len(parts) >= 6:  # Check if the line has at least 6 columns
        # Remove 'chr' from the chromosome column
        parts[1] = parts[1].replace('chr', '')
        # Limiting the third and fourth bit after ':' to 10 characters max
        split_chr = parts[1].split(':')
        if len(split_chr) >= 4:
            split_chr[2] = split_chr[2][:10]
            split_chr[3] = split_chr[3][:10]
            parts[1] = ':'.join(split_chr)
        return '\t'.join(parts)
    else:
        print("Warning: Skipping line with fewer columns than expected:", line.strip())
        return None

def main(input_file, output_file):
    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        header = next(f_in)  # Read the header line
        f_out.write(header)  # Write the header line to the output file
        for line in f_in:
            processed_line = process_line(line)
            if processed_line:
                f_out.write(processed_line + '\n')

    print("Processing complete. Output written to", output_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process whitespace-delimited file')
    parser.add_argument('-i', '--input', help='Input file path', required=True)
    parser.add_argument('-o', '--output', help='Output file path', required=True)
    args = parser.parse_args()
    
    main(args.input, args.output)
