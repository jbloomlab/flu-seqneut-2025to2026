# convert_csv_to_fa.py
# Example usage: python script.py input.csv output.fa`


import csv
import sys

if len(sys.argv) != 3:
    print("Usage: python script.py <input.csv> <output.fa>")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

with open(input_file, 'r') as csv_file, open(output_file, 'w') as fasta_file:
    reader = csv.DictReader(csv_file)
    seen = set()  # Track strain/sequence combinations we've already written
    
    for row in reader:
        # Clean the sequence: remove spaces, newlines, and other whitespace
        clean_sequence = ''.join(row['protein_sequence_HA_ectodomain'].split())
        strain = row['strain']
        
        # Create a unique identifier for this strain/sequence combination
        strain_seq_combo = (strain, clean_sequence)
        
        # Only write if we haven't seen this combination before
        if strain_seq_combo not in seen:
            seen.add(strain_seq_combo)
            fasta_file.write(f">{strain}\n")
            fasta_file.write(f"{clean_sequence}\n")

print(f"FASTA file '{output_file}' created successfully! Duplicates removed.")