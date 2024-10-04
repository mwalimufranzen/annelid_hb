from Bio import SeqIO
from Bio.Seq import Seq

# Input files
genome_file = "path_to_your_genome_file.fasta"  # Replace with the path to your genome file
shotgun_sequence_file = "path_to_shotgun_sequence.fasta"  # Replace with the path to your shotgun sequence

# Load the shotgun sequence from the file
def load_shotgun_sequence(file_path):
    with open(file_path, "r") as handle:
        shotgun_seq_record = SeqIO.read(handle, "fasta")
        return str(shotgun_seq_record.seq)

# Function to search for the shotgun sequence and its reverse complement in the genome
def search_in_genome(genome_file, shotgun_sequence):
    reverse_complement = str(Seq(shotgun_sequence).reverse_complement())
    
    for record in SeqIO.parse(genome_file, "fasta"):
        genome_seq = str(record.seq)
        
        # Search for the forward sequence
        if shotgun_sequence in genome_seq:
            start_pos = genome_seq.find(shotgun_sequence)
            end_pos = start_pos + len(shotgun_sequence)
            print(f"Found forward sequence in {record.id}: Start = {start_pos}, End = {end_pos}")
            return record.id, start_pos, end_pos, "forward"
        
        # Search for the reverse complement sequence
        if reverse_complement in genome_seq:
            start_pos = genome_seq.find(reverse_complement)
            end_pos = start_pos + len(reverse_complement)
            print(f"Found reverse complement in {record.id}: Start = {start_pos}, End = {end_pos}")
            return record.id, start_pos, end_pos, "reverse complement"

    print("Shotgun sequence not found in genome.")
    return None, None, None, None

# Main script
if __name__ == "__main__":
    shotgun_sequence = load_shotgun_sequence(shotgun_sequence_file)
    scaffold, start, end, orientation = search_in_genome(genome_file, shotgun_sequence)
    if scaffold:
        print(f"Shotgun sequence found in scaffold {scaffold} from position {start} to {end} ({orientation}).")
