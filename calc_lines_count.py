import gzip

# Initialize a counter for the total number of lines
total_lines = 0

# Loop through all the files
for i in range(1, 22):
    file_path = f"Colon_Sigmoid.v8.EUR.allpairs.chr{i}.txt.gz"
    
    # Open each file and count the number of lines
    with gzip.open(file_path, 'rt') as f:
        for _ in f:
            total_lines += 1

print(f"Total number of lines across all files: {total_lines}")
