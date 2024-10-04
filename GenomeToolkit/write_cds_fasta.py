import re

# Define input and output file paths
input_file_path = "alignment.txt"
output_file_path = ""

# Initialize variables
protein_id = ""
amino_acid_sequence = []

# Open the alignment file and process it
with open(input_file_path, 'r') as file:
    for line in file:
        # Look for the Protein ID line and extract the ID
        if line.startswith("Protein ID:"):
            protein_id_match = re.search(r"Protein ID: (\S+)", line)
            if protein_id_match:
                protein_id = protein_id_match.group(1)
        
        # Look for lines that start with 'Sbjct' and extract the amino acid sequence
        if line.startswith("Sbjct"):
            # Remove the numbers and spaces, keeping only the amino acid sequence
            sequence_match = re.search(r"Sbjct\s+\d+\s+([A-Z]+)\s+\d+", line)
            if sequence_match:
                amino_acid_sequence.append(sequence_match.group(1))

# Combine all the amino acid sequences into a single string
amino_acid_sequence_str = ''.join(amino_acid_sequence)

# Set the output file name based on the protein ID
output_file_path = f"{protein_id}_CDS.fasta"

# Write the amino acid sequence to the FASTA file with 60 characters per line
with open(output_file_path, 'w') as fasta_file:
    # Write the header
    fasta_file.write(f">{protein_id}_CDS.fasta\n")
    
    # Write the sequence in lines of 60 characters
    for i in range(0, len(amino_acid_sequence_str), 60):
        fasta_file.write(amino_acid_sequence_str[i:i+60] + "\n")

print(f"FASTA file written to {output_file_path}")
