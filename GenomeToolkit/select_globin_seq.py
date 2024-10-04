import re

def extract_globin_sequences_with_species(blast_output_file, fasta_output_file, min_length=120, max_length=180):
    # Define keywords to search for in the sequence descriptions
    keywords = ["Globin", "Myoglobin", "Hemoglobin", "Gb"]

    # Open the BLAST output file and read the content
    with open(blast_output_file, 'r') as blast_file:
        blast_data = blast_file.read()

    # Find headers that match the keywords and capture sequence IDs and species names
    # The regex pattern captures: >Sequence_ID Description text [Genus species]
    headers = re.findall(r'>([^\s]+)[^\n]*?(' + '|'.join(keywords) + r')[^\n]*?\[([A-Za-z]+ [a-z]+)\]', blast_data, re.IGNORECASE)

    if not headers:
        print("No matching sequences found with the specified keywords.")
        return

    # Find all Sbjct lines and capture the amino acid sequences
    sbjct_sequences = re.findall(r'Sbjct\s+[\d]+\s+([A-Za-z-]+)\s+[\d]+', blast_data)

    # Create a dictionary to map headers (ID + species) to their sequences
    header_to_sequence = {}
    header_indices = [m.start() for m in re.finditer(r'>', blast_data)]

    # Iterate over headers to capture corresponding Sbjct sequences
    for idx, (header_id, _, species) in enumerate(headers):
        # Get the index range in the BLAST data for the current header
        start_idx = header_indices[idx]
        end_idx = header_indices[idx + 1] if idx + 1 < len(header_indices) else len(blast_data)

        # Extract all Sbjct sequences for this particular header
        current_section = blast_data[start_idx:end_idx]
        sbjct_lines = re.findall(r'Sbjct\s+[\d]+\s+([A-Za-z-]+)\s+[\d]+', current_section)

        # Concatenate all lines into a single sequence and remove gaps
        full_sequence = ''.join(sbjct_lines).replace("-", "")
        
        # Filter based on sequence length (between min_length and max_length)
        if min_length <= len(full_sequence) <= max_length:
            # Construct the FASTA header with sequence ID and species name
            header_to_sequence[f"{header_id} [{species}]"] = full_sequence

    # Write the results to a FASTA file
    with open(fasta_output_file, 'w') as fasta_file:
        for header, sequence in header_to_sequence.items():
            # Write the header in FASTA format
            fasta_file.write(f">{header}\n")
            # Write the sequence, breaking into 60 characters per line
            for i in range(0, len(sequence), 60):
                fasta_file.write(sequence[i:i+60] + '\n')

    print(f"Extracted {len(header_to_sequence)} filtered globin sequences to {fasta_output_file}.")

# File names (input and output)
blast_output = "blastp_outfmt0.txt"  # Change to your actual BLAST output file
fasta_output = "extracted_blastp_seq.fasta"

# Run the function to extract filtered globin sequences with species information
extract_globin_sequences_with_species(blast_output, fasta_output, min_length=90, max_length=180)
