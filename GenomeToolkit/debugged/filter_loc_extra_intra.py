from Bio import SeqIO

def filter_and_categorize_fasta(input_fasta, min_length=110, max_length=180):
    """
    Filters sequences in a FASTA file by length and categorizes them based on whether they are 
    extracellular or non-extracellular.

    Args:
        input_fasta (str): Path to the input multi-FASTA file.
        min_length (int): Minimum length of sequences to retain.
        max_length (int): Maximum length of sequences to retain.
        
    Output:
        Writes two output FASTA files: one for extracellular globins and one for non-extracellular globins.
    """
    extracellular_file = "extracellular_globins.fasta"
    non_extracellular_file = "non_extracellular_globins.fasta"

    # Create file handles for the two categories
    extracellular_records = []
    non_extracellular_records = []

    # Read through the input FASTA and classify sequences
    for record in SeqIO.parse(input_fasta, "fasta"):
        sequence_length = len(record.seq)
        
        # Filter based on length first
        if min_length <= sequence_length <= max_length:
            # Check if the description contains "extracellular"
            if any(keyword in record.description.lower() for keyword in ["extracellular", "secreted", "plasma"]):
                extracellular_records.append(record)
            else:
                non_extracellular_records.append(record)

#            if "extracellular" in record.description.lower():
#               extracellular_records.append(record)
#            else:
#                non_extracellular_records.append(record)

    # Write sequences to the respective files
    SeqIO.write(extracellular_records, extracellular_file, "fasta")
    SeqIO.write(non_extracellular_records, non_extracellular_file, "fasta")

    # Output the results
    print(f"Filtered and categorized sequences have been saved to:\n"
          f"  1. '{extracellular_file}' (Extracellular Globins): {len(extracellular_records)} sequences\n"
          f"  2. '{non_extracellular_file}' (Non-Extracellular Globins): {len(non_extracellular_records)} sequences")


# Specify the input file path
# Specify non-redundant globins filtered using cd-hit
input_file = "annelid_globins_nr.fasta"  # Replace with your actual input file name

# Run the filtering and categorization function
filter_and_categorize_fasta(input_file)
