from Bio import SeqIO
from collections import defaultdict

def parse_fasta_files(fasta_files):
    """
    Parse multiple FASTA files and extract species information, sequence counts, and average sequence length.

    Args:
        fasta_files (list): List of FASTA file paths to be read.

    Returns:
        dict: A dictionary where keys are species names and values are tuples (globin_count, total_length).
    """
    species_data = defaultdict(lambda: [0, 0])  # Stores species_name -> [globin_count, total_length]

    # Process each FASTA file and accumulate species data
    for fasta_file in fasta_files:
        for record in SeqIO.parse(fasta_file, "fasta"):
            # Extract the species name, which is found in brackets [Genus species]
            if "[" in record.description and "]" in record.description:
                species_name = record.description.split("[")[1].split("]")[0].strip()
                
                # Update the count and total length for each species
                species_data[species_name][0] += 1  # Increment globin count
                species_data[species_name][1] += len(record.seq)  # Accumulate total sequence length

    return species_data

def write_species_summary(species_data, output_file):
    """
    Write the species summary to a file.

    Args:
        species_data (dict): Dictionary with species names as keys and (globin_count, total_length) as values.
        output_file (str): Path to the output file.
    """
    with open(output_file, "w") as f:
        # Write header
        f.write("Species\tNumber_of_Globins\tAverage_Length\n")

        # Write each species and its corresponding data
        for species, (globin_count, total_length) in sorted(species_data.items()):
            avg_length = total_length / globin_count if globin_count > 0 else 0
            f.write(f"{species}\t{globin_count}\t{avg_length:.2f}\n")

    print(f"Species summary written to '{output_file}'.")

# Specify the input FASTA files
input_fasta_files = ["extracellular_globins.fasta", "non_extracellular_globins.fasta"]

# Specify the output file for the species summary
output_species_summary = "species_globin_summary.txt"

# Parse the input FASTA files and gather species data
species_info = parse_fasta_files(input_fasta_files)

# Write the species summary to the output file
write_species_summary(species_info, output_species_summary)
