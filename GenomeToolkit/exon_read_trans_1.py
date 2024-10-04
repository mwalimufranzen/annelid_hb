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
                orientation = "sense"
            elif sbjct_start == end:
                orientation = "antisense"
            else:
                sys.stderr.write(f"Unable to determine strand orientation for range {start}-{end}.\n")
                sys.exit(1)
            ranges.append((start, end, orientation))
        else:
            ranges.append((start, end, 'unknown'))  # Allow to proceed even without orientation

    return protein_id, sequence_id, ranges

def concatenate_exons(exon_files):
    concatenated_seq = ""
    for exon_file in exon_files:
        exon_seq = read_fasta_file(exon_file)
        if exon_seq:  # Only concatenate if exon_seq is not None
            concatenated_seq += str(exon_seq)
    return concatenated_seq

def translate_sequence(nucleotide_seq, orientation):
    """Translate the sequence based on the orientation (sense or antisense)"""
    orf = None
    if orientation == "sense" and nucleotide_seq[:3] == "ATG":
        orf = nucleotide_seq  # First reading frame if start codon is ATG
    elif orientation == "antisense" and nucleotide_seq[-3:] == "CAT":
        orf = nucleotide_seq.reverse_complement()  # Fourth reading frame if reverse complement
    else:
        sys.stderr.write("Error: No valid start codon (ATG or CAT) found.\n")
        sys.exit(1)

    # Translate the ORF sequence
    amino_acid_seq = str(Seq(orf).translate(to_stop=True))
    return amino_acid_seq

def write_fasta(seq, file_name, description):
    seq_record = SeqRecord(Seq(seq), id="coding_region", description=description)
    with open(file_name, 'w') as fasta_out:
        SeqIO.write(seq_record, fasta_out, 'fasta')

def main():
    # Parse the alignment file to get information about the exons and orientation
    alignment_file = "alignment.txt"
    protein_id, sequence_id, exon_positions = parse_alignment_file(alignment_file)

    # Input exon files (change as per available files)
    exon_files = ['exon1.fasta', 'exon2.fasta', 'exon3.fasta', 'exon4.fasta']
    exon_files = [file for file in exon_files if os.path.exists(file)]  # Filter out non-existent files

    # Step 1: Concatenate exons
    concatenated_seq = concatenate_exons(exon_files)
    if concatenated_seq:
        write_fasta(concatenated_seq, f"{protein_id}_coding.fasta", f"{sequence_id} ORF for {protein_id}")

        # Step 2: Translate the coding sequence based on orientation
        if exon_positions:
            orientation = exon_positions[0][2]  # Assuming all exons share the same orientation
            print(f"Orientation: {orientation}")

            # Translate sequence and write output
            translated_seq = translate_sequence(concatenated_seq, orientation)
            write_fasta(translated_seq, f"{protein_id}_trans.fasta", f"{sequence_id} {orientation}ORF for {protein_id}")
        else:
            sys.stderr.write("Error: No exon positions found in the alignment file.\n")
    else:
        print("No exons found to concatenate.")

if __name__ == "__main__":
    main()
