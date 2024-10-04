import re

def filter_annelid_sequences(input_fasta, output_fasta):
    # Define a list of annelid genera or species
    annelid_genera = [
        "Capitella", "Alitta", "Hirudo", "Lumbricus", "Nereis", "Branchipolynoe", 
        "Owenia", "Spio", "Arenicola", "Platynereis", "Terebellum", "Serpula"
    ]
    annelid_species = [
        "Capitella teleta", "Alitta virens", "Hirudo medicinalis", "Lumbricus terrestris",
        "Nereis diversicolor", "Branchipolynoe seepensis", "Owenia fusiformis"
    ]

    # Open the input FASTA file and read the content
    with open(input_fasta, 'r') as file:
        fasta_data = file.read()

    # Split the data into entries (each entry starts with a '>' line)
    entries = fasta_data.split('>')[1:]  # Split and ignore the empty first split before the first '>'
    
    # Store filtered annelid sequences
    annelid_sequences = []

    print("Starting to process FASTA entries...")

    # Process each FASTA entry
    for entry in entries:
        # Extract the header line and the sequence
        lines = entry.strip().split('\n')
        header = lines[0]
        sequence = ''.join(lines[1:])  # Join all sequence lines into a single string

        # Extract genus and species information from the header using regex
        match = re.search(r'\[([A-Za-z]+ [a-z]+)\]', header)
        if match:
            species_info = match.group(1)  # Get the species name like 'Alitta virens'
            genus = species_info.split()[0]  # Extract the genus from 'Alitta virens'

            print(f"Processing sequence: {header}")  # Print the current header being processed
            print(f"Length of sequence: {len(sequence)} amino acids")  # Print the length of the sequence

            # Check if the genus or species is in the annelid list
            if genus in annelid_genera or species_info in annelid_species:
                # Print confirmation of annelid match
                print(f"Annelid match found: {genus} [{species_info}]")

                # Add the matching entry to the filtered list
                formatted_sequence = "\n".join([sequence[i:i+60] for i in range(0, len(sequence), 60)])
                annelid_sequences.append(f">{header}\n{formatted_sequence}\n")
            else:
                print(f"Not an annelid: {genus} [{species_info}]")

    # Write the results to a new FASTA file if any annelid sequences were found
    if annelid_sequences:
        with open(output_fasta, 'w') as output_file:
            output_file.writelines(annelid_sequences)
        print(f"Filtered {len(annelid_sequences)} annelid sequences to {output_fasta}.")
    else:
        print("No annelid sequences found in the input file.")

# Input and output FASTA files
input_fasta_file = "filtered_globin_sequences_with_species.fasta"
output_fasta_file = "annelid_globin_sequences.fasta"

# Run the filtering function with print statements for debugging
filter_annelid_sequences(input_fasta_file, output_fasta_file)
