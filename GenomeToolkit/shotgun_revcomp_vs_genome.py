from Bio import SeqIO
from Bio.Seq import Seq

# Input files
genome_file = "path_to_your_genome_file.fasta"  # Replace with the path to your genome file
shotgun_sequence_file = "path_to_shotgun_sequence.fasta"  # Replace with the path to your shotgun sequence
output_file = "genome_match_region.fasta"  # Output file to save the matched region

# Load the shotgun sequence from the file
def load_shotgun_sequence(file_path):
    with open(file_path, "r") as handle:
        shotgun_seq_record = SeqIO.read(handle, "fasta")
        return str(shotgun_seq_record.seq)

# Function to search for the shotgun sequence and its reverse complement in the genome
def search_in_genome(genome_file, shotgun_sequence, flanking_size=4000):
    reverse_complement = str(Seq(shotgun_sequence).reverse_complement())
    
    for record in SeqIO.parse(genome_file, "fasta"):
        genome_seq = str(record.seq)
        
        # Search for the forward sequence
        if shotgun_sequence in genome_seq:
            start_pos = genome_seq.find(shotgun_sequence)
            end_pos = start_pos + len(shotgun_sequence)
            print(f"Found forward sequence in {record.id}: Start = {start_pos}, End = {end_pos}")
            return record, start_pos, end_pos, "forward", flanking_size
        
        # Search for the reverse complement sequence
        if reverse_complement in genome_seq:
            start_pos = genome_seq.find(reverse_complement)
            end_pos = start_pos + len(reverse_complement)
            print(f"Found reverse complement in {record.id}: Start = {start_pos}, End = {end_pos}")
            return record, start_pos, end_pos, "reverse complement", flanking_size

    print("Shotgun sequence not found in genome.")
    return None, None, None, None, None

# Function to extract the matched region and flanking nucleotides
def extract_flanking_region(record, start_pos, end_pos, flanking_size, orientation, output_file):
    genome_seq = record.seq
    
    # Extract flanking region, with boundary checks
    flanking_start = max(0, start_pos - flanking_size)
    flanking_end = min(len(genome_seq), end_pos + flanking_size)
    
    match_region = genome_seq[flanking_start:flanking_end]
    
    # If it's the reverse complement, flip the orientation
    if orientation == "reverse complement":
        match_region = match_region.reverse_complement()

    # Write the flanking region to a FASTA file
    with open(output_file, "w") as output_handle:
        SeqIO.write(SeqIO.SeqRecord(match_region, id=f"{record.id}_region", description=f"flanking_{flanking_size}_nts"), output_handle, "fasta")
    
    print(f"Match region with flanking nucleotides written to {output_file}")

# Main script
if __name__ == "__main__":
    shotgun_sequence = load_shotgun_sequence(shotgun_sequence_file)
    record, start, end, orientation, flanking_size = search_in_genome(genome_file, shotgun_sequence)
    if record:
        extract_flanking_region(record, start, end, flanking_size, orientation, output_file)
        print(f"Shotgun sequence found in {record.id} (scaffold/chromosome) from position {start} to {end} ({orientation}).")
