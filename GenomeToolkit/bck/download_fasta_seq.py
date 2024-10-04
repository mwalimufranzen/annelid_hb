from Bio import Entrez
import re
import sys

def extract_sequence_id(alignment_file):
    with open(alignment_file, 'r') as file:
        content = file.read()
        # Extract sequence ID from the alignment file
        match_id = re.search(r'Sequence ID:\s*(\S+)', content)
        if match_id:
            return match_id.group(1).strip()
        else:
            sys.stderr.write("No sequence ID found in the alignment file.\n")
            sys.exit(1)

def fetch_and_save_sequence(sequence_id, output_file):
    Entrez.email = 'franzen@ncsu.edu'  # Replace with your email
    try:
        # Fetch the sequence from NCBI
        sys.stderr.write(f"Fetching sequence for ID: {sequence_id}\n")
        handle = Entrez.efetch(db="nuccore", id=sequence_id, rettype="fasta", retmode="text")
        sequence_data = handle.read()
        
        # Write the sequence data to a file
        with open(output_file, 'w') as file:
            file.write(sequence_data)
        
        sys.stderr.write(f"FASTA file written to: {output_file}\n")
    except Exception as e:
        sys.stderr.write(f"Error fetching or writing sequence: {e}\n")
        sys.exit(1)

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Fetch a FASTA sequence from NCBI based on an alignment file.")
    parser.add_argument('alignment_file', help="File containing alignment information including the sequence ID")
    parser.add_argument('output_file', help="Output file to save the FASTA sequence")
    
    args = parser.parse_args()
    
    sequence_id = extract_sequence_id(args.alignment_file)
    fetch_and_save_sequence(sequence_id, args.output_file)

