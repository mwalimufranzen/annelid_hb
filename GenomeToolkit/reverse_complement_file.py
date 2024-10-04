from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def reverse_complement(sequence):
    """Return the reverse complement of the given sequence."""
    seq = Seq(sequence)
    return str(seq.reverse_complement())

def format_fasta(sequence, line_length=70):
    """Format the sequence in FASTA format with lines of specified length."""
    return '\n'.join(sequence[i:i+line_length] for i in range(0, len(sequence), line_length))

def process_fasta(input_file, output_file):
    """Read a FASTA file, compute the reverse complement of the sequences, and write to a new FASTA file with formatted headers and sequences."""
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        for record in SeqIO.parse(infile, "fasta"):
            # Compute reverse complement
            rev_comp_seq = reverse_complement(str(record.seq))
            # Reverse the order of elements in the header
            header_parts = record.id.split('|')
            reversed_header = '|'.join(reversed(header_parts))
            # Format the sequence in lines of 70 characters
            formatted_seq = format_fasta(rev_comp_seq, line_length=70)
            # Write to output file
            outfile.write(f">{reversed_header}\n{formatted_seq}\n")

if __name__ == "__main__":
    # Define input and output file names
    input_fasta = "seq2_rc.fasta"
    output_fasta = "seq2.fasta"
    
    # Process the FASTA file
    process_fasta(input_fasta, output_fasta)
