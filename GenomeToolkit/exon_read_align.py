from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import PairwiseAligner
import re
import os

def read_fasta_file(file_name):
    """Reads a single sequence from a FASTA file."""
    if os.path.exists(file_name):
        with open(file_name, 'r') as fasta_file:
            for record in SeqIO.parse(fasta_file, 'fasta'):
                return record.seq
    return None

def parse_alignment_file(alignment_file):
    """Parse the alignment.txt file to extract Protein ID, Sequence ID, exon positions, and amino acid sequence."""
    protein_id = None
    amino_acid_sequence = []
    orientation = None
    range_start, range_end = None, None

    with open(alignment_file, 'r') as file:
        for line in file:
            # Extract the protein ID
            if line.startswith("Protein ID:"):
                protein_id_match = re.search(r"Protein ID: (\S+)", line)
                if protein_id_match:
                    protein_id = protein_id_match.group(1)
            
            # Extract the range from the first occurrence of 'Range'
            if line.startswith("Range"):
                range_match = re.search(r"Range \d+: (\d+) to (\d+)", line)
                if range_match:
                    range_start, range_end = map(int, range_match.groups())
            
            # Detect orientation by comparing Sbjct number with range values
            if line.startswith("Sbjct") and range_start is not None and range_end is not None:
                sbjct_match = re.search(r"Sbjct\s+(\d+)\s+[A-Z]+\s+(\d+)", line)
                if sbjct_match:
                    sbjct_start = int(sbjct_match.group(1))
                    # Compare Sbjct start with range_start and range_end to determine orientation
                    if sbjct_start == range_start:
                        #print(f"Detected sense orientation")
                        orientation = "+"
                    elif sbjct_start == range_end:
                        #print(f"Detected antisense orientation")
                        orientation = "-"
            
            # Extract Sbjct amino acid sequences
            if line.startswith("Sbjct"):
                sequence_match = re.search(r"Sbjct\s+\d+\s+([A-Z]+)\s+\d+", line)
                if sequence_match:
                    amino_acid_sequence.append(sequence_match.group(1))
    
    # Combine all the amino acid sequences
    amino_acid_sequence_str = ''.join(amino_acid_sequence)
    
    return protein_id, amino_acid_sequence_str, orientation

def concatenate_exons(exon_files, orientation):
    """Concatenate exons in the correct order based on orientation."""
    concatenated_seq = ""
    
    if orientation == "-":
        exon_files = exon_files[::-1]  # Reverse order for antisense

    for exon_file in exon_files:
        exon_seq = read_fasta_file(exon_file)
        if exon_seq:
            concatenated_seq += str(exon_seq)
    
    return concatenated_seq

def translate_sequence(nucleotide_seq, orientation):
    """Translate the nucleotide sequence to amino acids."""
    nucleotide_seq = Seq(nucleotide_seq)
    
    if orientation == "-":
        print(f"Taking the reverse complement of the sequence for translation.")
        nucleotide_seq = nucleotide_seq.reverse_complement()  # Reverse complement for antisense
    
    return str(nucleotide_seq.translate())

def write_fasta(seq, file_name, description):
    """Write a sequence to a FASTA file."""
    seq_record = SeqRecord(Seq(seq), description=description, id="coding_region")
    with open(file_name, 'w') as fasta_out:
        SeqIO.write(seq_record, fasta_out, 'fasta')

def run_pairwise_alignment(seq1, seq2, output_file="comparison.txt"):
    """Perform pairwise alignment and write the result to a file."""
    aligner = PairwiseAligner()
    alignment = aligner.align(seq1, seq2)
    with open(output_file, "w") as output:
        output.write("Pairwise Alignment:\n")
        output.write(str(alignment[0]))  # Write the first alignment only

def main():
    alignment_file = "alignment.txt"
    
    # Step 1: Parse the alignment file to get Protein ID, CDS, and amino acid sequence
    protein_id, cds_aa_sequence, orientation = parse_alignment_file(alignment_file)
    
    # Step 2: Write the CDS FASTA file (from alignment)
    cds_fasta_file = f"{protein_id}_CDS.fasta"
    write_fasta(cds_aa_sequence, cds_fasta_file, f"CDS from alignment for {protein_id}")
    
    # Step 3: Concatenate exons from exon FASTA files
    print(f"Concatenating 3 exon files. Alter code for other number of files.")
    exon_files = ['exon1.fasta', 'exon2.fasta', 'exon3.fasta']  # Replace with actual exon files
    concatenated_exon_seq = concatenate_exons(exon_files, orientation)
    
    # Step 4: Write the concatenated exon sequence to a FASTA file
    coding_fasta_file = f"{protein_id}_coding.fasta"
    write_fasta(concatenated_exon_seq, coding_fasta_file, f"Sequence for {protein_id}")
    
    # Step 5: Translate the concatenated exon sequence
    translated_seq = translate_sequence(concatenated_exon_seq, orientation)
    
    # Step 6: Write the translated sequence to a FASTA file
    trans_fasta_file = f"{protein_id}_trans.fasta"
    write_fasta(translated_seq, trans_fasta_file, f"Translated sequence for {protein_id}")
    
    # Step 7: Compare the translated sequence with the CDS from the alignment file
    alignment_output_file = f"{protein_id}_align.txt"
    print(f"Aligning target {protein_id}_CDS with query {protein_id}_trans.")
    run_pairwise_alignment(cds_aa_sequence, translated_seq, alignment_output_file)

if __name__ == "__main__":
    main()
