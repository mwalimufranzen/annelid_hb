from Bio import SeqIO

def filter_fasta_by_length(input_fasta, output_fasta, min_length=110, max_length=180):
    """
    Filters sequences in a FASTA file by length and writes the filtered sequences to a new file.
    
    Args:
        input_fasta (str): Path to the input FASTA file.
        output_fasta (str): Path to the output FASTA file.
        min_length (int): Minimum length of sequences to retain.
        max_length (int): Maximum length of sequences to retain.
    """
    # Read the input FASTA file and filter based on sequence length
    filtered_sequences = (record for record in SeqIO.parse(input_fasta, "fasta")
                          if min_length <= len(record.seq) <= max_length)

    # Write the filtered sequences to the output file
    count = SeqIO.write(filtered_sequences, output_fasta, "fasta")
    
    # Print out how many sequences were written
    print(f"Filtered {count} sequences between {min_length} and {max_length} amino acids long and saved to '{output_fasta}'.")

# Specify input and output file paths
input_file = "input_sequences.fasta"  # Replace with your actual input file
output_file = "filtered_sequences.fasta"  # Replace with your desired output file

# Run the filtering function
filter_fasta_by_length(input_file, output_file)
