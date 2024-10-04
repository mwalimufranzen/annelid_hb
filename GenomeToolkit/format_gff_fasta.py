def calculate_sequence_lengths(fasta_file):
    """
    Reads the sequences from the FASTA file and calculates their lengths.
    Returns a dictionary with sequence ID as the key and the length of the sequence as the value.
    """
    sequence_lengths = {}
    sequence_id = None
    sequence = []
    
    print(f"Reading FASTA from {fasta_file}")
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()

            if line.startswith(">"):
                if sequence_id is not None:
                    # Store the previous sequence length
                    sequence_lengths[sequence_id] = len(''.join(sequence))
                    print(f"Sequence ID: {sequence_id}, Length: {sequence_lengths[sequence_id]}")
                # Start new sequence
                sequence_id = line[1:].split()[0]
                print(f"Found sequence: {sequence_id}")
                sequence = []
            else:
                # Append sequence data
                sequence.append(line)
        
        # Store the last sequence length
        if sequence_id is not None:
            sequence_lengths[sequence_id] = len(''.join(sequence))
            print(f"Sequence ID: {sequence_id}, Length: {sequence_lengths[sequence_id]}")

    return sequence_lengths

def add_sequence_region_statements(gff_file, fasta_file, output_file):
    """
    Adds ##sequence-region statements to the GFF file based on the lengths of the sequences in the FASTA file.
    """
    sequence_lengths = calculate_sequence_lengths(fasta_file)
    seen_sequences = set()
    
    with open(gff_file, 'r') as gff, open(output_file, 'w') as output:
        print(f"Processing GFF file {gff_file}")

        for line in gff:
            line = line.strip()
            
            # Process the GFF section
            fields = line.split()
            if len(fields) > 0 and fields[0] in sequence_lengths:
                seq_id = fields[0]
                if seq_id not in seen_sequences:
                    # Add a new ##sequence-region line when sequence ID changes
                    output.write(f"##sequence-region {seq_id} 1 {sequence_lengths[seq_id]}\n")
                    print(f"Writing ##sequence-region for {seq_id}")
                    seen_sequences.add(seq_id)
            
            # Write the current GFF line
            output.write(line + '\n')

# Usage
gff_file = "Annelid_Hb.gff"  # Path to your input GFF file
fasta_file = "Annelid_Hb_seq.fasta"  # Path to your FASTA file
output_file = "Annelid_Hb_with_regions.gff"  # Output file

add_sequence_region_statements(gff_file, fasta_file, output_file)
