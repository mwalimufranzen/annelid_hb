from Bio import SeqIO
from collections import defaultdict

def parse_fasta_files(fasta_files):
    """
    Parse multiple FASTA files and extract genus and species information, sequence counts, 
    total sequence length, and amino acid percentages.

    Args:
        fasta_files (list): List of FASTA file paths to be read.

    Returns:
        dict: A dictionary where keys are 'Genus species' and values are tuples
              (globin_count, total_length, total_aa_percentages).
    """
    # Stores 'Genus species' -> [globin_count, total_length, aa_sums]
    species_data = defaultdict(lambda: [0, 0, {'M': 0, 'Y': 0, 'H': 0, 'C': 0, 'W': 0}])

    # Process each FASTA file and accumulate species data
    for fasta_file in fasta_files:
        for record in SeqIO.parse(fasta_file, "fasta"):
            if "OS=" in record.description:
                # Extract the species name from the description using the 'OS=' pattern
                species_name = " ".join(record.description.split("OS=")[1].split(" ")[:2]).strip()

                # Update the count and total length for each species
                sequence = str(record.seq)
                species_data[species_name][0] += 1  # Increment globin count
                species_data[species_name][1] += len(sequence)  # Accumulate total sequence length

                # Calculate amino acid percentages
                aa_percentages = calculate_amino_acid_percentages(sequence)
                for aa, percent in aa_percentages.items():
                    species_data[species_name][2][aa] += percent  # Sum up the percentages

    return species_data

def calculate_amino_acid_percentages(seq):
    """
    Calculate the percentages of Methionine (M), Tyrosine (Y), Histidine (H),
    Cysteine (C), and Tryptophan (W) in a given sequence.

    Args:
        seq (str): The amino acid sequence.

    Returns:
        dict: A dictionary of amino acid percentages.
    """
    total_length = len(seq)
    counts = {
        'M': seq.count('M'),
        'Y': seq.count('Y'),
        'H': seq.count('H'),
        'C': seq.count('C'),
        'W': seq.count('W')
    }
    percentages = {aa: (count / total_length) * 100 for aa, count in counts.items()}
    return percentages

def write_species_summary(species_data, output_file):
    """
    Write the species summary with amino acid percentages to a file.

    Args:
        species_data (dict): Dictionary with species names as keys and (globin_count, total_length, aa_sums) as values.
        output_file (str): Path to the output file.
    """
    with open(output_file, "w") as f:
        # Write header
        f.write("Genus_Species\tNumber_of_Globins\tAverage_Length\tMet(%)\tTyr(%)\tHis(%)\tCys(%)\tTrp(%)\n")

        # Write each species and its corresponding data
        for species, (globin_count, total_length, aa_sums) in sorted(species_data.items()):
            avg_length = total_length / globin_count if globin_count > 0 else 0
            avg_aa_percentages = {aa: aa_sums[aa] / globin_count if globin_count > 0 else 0 for aa in aa_sums}
            f.write(f"{species}\t{globin_count}\t{avg_length:.2f}\t"
                    f"{avg_aa_percentages['M']:.2f}\t{avg_aa_percentages['Y']:.2f}\t"
                    f"{avg_aa_percentages['H']:.2f}\t{avg_aa_percentages['C']:.2f}\t"
                    f"{avg_aa_percentages['W']:.2f}\n")

    print(f"Species summary with amino acid percentages written to '{output_file}'.")

# Specify the input FASTA files
input_fasta_files = ["extracellular_globins.fasta", "non_extracellular_globins.fasta"]

# Specify the output file for the species summary
output_species_summary = "species_globin_summary.txt"

# Parse the input FASTA files and gather species data
species_info = parse_fasta_files(input_fasta_files)

# Write the species summary to the output file
write_species_summary(species_info, output_species_summary)
