from Bio import SeqIO
from Bio.Align import PairwiseAligner

def read_fasta_sequence(fasta_file):
    """Reads the first sequence from a FASTA file."""
    with open(fasta_file, "r") as file:
        for record in SeqIO.parse(file, "fasta"):
            return str(record.seq)

def run_pairwise_alignment(seq1, seq2):
    """Performs pairwise alignment using Biopython's PairwiseAligner."""
    aligner = PairwiseAligner()
    aligner.mode = 'global'  # Use global alignment, can be set to 'local' if needed
    alignment = aligner.align(seq1, seq2)
    return alignment

def compare_sequences_with_alignment(seq1, seq2, output_file="comparison.txt"):
    """Performs pairwise alignment on two sequences and writes the first alignment to a file."""
    
    # Run Biopython's PairwiseAligner
    alignment = run_pairwise_alignment(seq1, seq2)
    
    # Write only the first alignment result to comparison.txt
    with open(output_file, "w") as output:
        output.write("Pairwise Alignment:\n")
        output.write(str(alignment[0]))  # Write the first alignment only
        output.write("\n\n")
    
    print(f"First alignment successfully written to {output_file}")

# Read sequences from the provided FASTA files
seq1 = read_fasta_sequence("seq1.fasta")  # query
seq2 = read_fasta_sequence("seq2.fasta") # target

# Compare the sequences and write the first alignment output to comparison.txt
compare_sequences_with_alignment(seq1, seq2)
