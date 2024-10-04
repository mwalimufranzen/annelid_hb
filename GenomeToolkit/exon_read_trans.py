from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import re
import sys

def read_fasta_file(file_name):
    if os.path.exists(file_name):  # Check if the file exists
        with open(file_name, 'r') as fasta_file:
            for record in SeqIO.parse(fasta_file, 'fasta'):
                return record.seq  # Assume single sequence per file
    else:
        return None  # Return None if the file does not exist


def parse_alignment_file(alignment_file):
    """Parse the alignment.txt file to extract relevant information."""
    exon_positions = []
    sequence_id = None
    protein_id = None
    content = None
    orientation = None

    # Read the alignment file content
    with open(alignment_file, 'r') as f:
        content = f.read()

    # Extract Protein ID
    match = re.search(r'Protein ID: (\w+)', content)
    if match:
        protein_id = match.group(1)
        print(f"Protein ID found: {protein_id}")

    # Extract Sequence ID (without "Length")
    match = re.search(r'Sequence ID: (\w+\.\d+)', content)
    if match:
        sequence_id = match.group(1)
        print(f"Sequence ID found: {sequence_id}")

    # Extract ranges and determine strand orientation based on the first occurrence of Sbjct line
    ranges = []
    for match in re.finditer(r'Range \d+: (\d+) to (\d+)', content):
        start, end = map(int, match.groups())
        sbjct_match = re.search(r'Sbjct\s+(\d+)\s+[^\n]+\s+(\d+)', content[match.end():])
        if sbjct_match:
            sbjct_start, sbjct_end = map(int, sbjct_match.groups())
            if sbjct_start == start:
                orientation = "+"  # Sense strand
            elif sbjct_start == end:
                orientation = "-"  # Antisense strand
            else:
                sys.stderr.write(f"Unable to determine strand orientation for range {start}-{end}.\n")
                sys.exit(1)
            ranges.append((start, end, orientation))
            sys.stdout.write(f"Range {len(ranges)}: start {start}, end {end}, {orientation}\n")
        else:
            sys.stderr.write(f"No Sbjct line found for range {start}-{end}. Trying to proceed.\n")
            ranges.append((start, end, 'unknown'))  # Allow it to proceed even without orientation

    # Store exon positions
    for idx, (start, end, orientation) in enumerate(ranges, 1):
        exon_name = f"exon{idx}"
        exon_positions.append((exon_name, start, end, orientation))

    return protein_id, sequence_id, exon_positions, orientation

def concatenate_exons(exon_files, orientation):
    concatenated_seq = ""
    
    # If orientation is antisense, reverse the exon order
    if orientation == "-":
        exon_files = exon_files[::-1]  # Reverse the order for antisense

    # Concatenate the exons in the correct order
    for exon_file in exon_files:
        exon_seq = read_fasta_file(exon_file)
        if exon_seq:  # Only concatenate if exon_seq is not None
            concatenated_seq += str(exon_seq)
    
    return concatenated_seq

def translate_sequence(nucleotide_seq, orientation):
    """Translate the sequence based on the orientation (sense or antisense)"""
    nucleotide_seq = Seq(nucleotide_seq)  # Ensure nucleotide_seq is a Seq object

    # Normalize orientation to '+' or '-'
    if orientation == "sense":
        orientation = "+"
    elif orientation == "antisense":
        orientation = "-"

    # Check orientation and handle translation
    if orientation == "+":
        orf = nucleotide_seq  # Use the sequence as is for sense strand
    elif orientation == "-":
        orf = nucleotide_seq.reverse_complement()  # Use the reverse complement for antisense strand
    else:
        raise ValueError("Invalid orientation: expected '+' or '-'")
    
    # Translate the ORF sequence, allowing '*' for stop codons
    amino_acid_seq = orf.translate()  # Removed 'to_stop=True' to allow '*' for stop codons
    
    return str(amino_acid_seq)

def write_fasta(seq, file_name, description):
    seq_record = SeqRecord(Seq(seq), description=description, id="coding_region")
    with open(file_name, 'w') as fasta_out:
        SeqIO.write(seq_record, fasta_out, 'fasta')

def main():
    alignment_file = "alignment.txt"
    
    # Parse the alignment file to get Protein ID, Sequence ID, exon positions, and orientation
    protein_id, sequence_id, exon_positions, orientation = parse_alignment_file(alignment_file)

    # Input exon files (ensure the correct exons are included)
    exon_files = ['exon1.fasta', 'exon2.fasta', 'exon3.fasta', 'exon4.fasta']
    exon_files = [file for file in exon_files if os.path.exists(file)]  # Filter out non-existent files

    if exon_positions:
        print(f"Orientation: {orientation}")

        # Step 1: Concatenate exons based on orientation
        concatenated_seq = concatenate_exons(exon_files, orientation)
        if concatenated_seq:
            write_fasta(concatenated_seq, f"{protein_id}_coding.fasta", f"{sequence_id} ORF for {protein_id}")

            # Step 2: Translate the coding sequence based on orientation
            translated_seq = translate_sequence(concatenated_seq, orientation)
            write_fasta(translated_seq, f"{protein_id}_trans.fasta", f"{sequence_id} {orientation} ORF for {protein_id}")
        else:
            print("No exons found to concatenate.")
    else:
        sys.stderr.write("Error: No exon positions found in the alignment file.\n")

if __name__ == "__main__":
    main()
