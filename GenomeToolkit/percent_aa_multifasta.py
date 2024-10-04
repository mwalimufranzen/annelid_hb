from Bio import SeqIO

def count_amino_acids(seq):
    """Count the specific amino acids in a sequence and calculate their percentages."""
    total_length = len(seq)
    counts = {
        'M': seq.count('M'),
        'Y': seq.count('Y'),
        'H': seq.count('H'),
        'C': seq.count('C'),
        'W': seq.count('W')
    }
    percentages = {aa: (count / total_length) * 100 for aa, count in counts.items()}
    return counts, percentages, total_length

def process_fasta(input_fasta, output_file):
    """Process the FASTA file and write the amino acid counts and percentages to the output file."""
    with open(output_file, 'w') as out:
        for record in SeqIO.parse(input_fasta, "fasta"):
            # Extract the accession number after 'tr|'
            accession = record.id.split('|')[1] if '|' in record.id else record.id
            # Count the amino acids and calculate percentages
            counts, percentages, total_length = count_amino_acids(str(record.seq))
            # Format the output line
            output_line = (
                f"{accession}\tLength: {total_length}\t"
                f"M: {counts['M']} ({percentages['M']:.2f}%)\t"
                f"Y: {counts['Y']} ({percentages['Y']:.2f}%)\t"
                f"H: {counts['H']} ({percentages['H']:.2f}%)\t"
                f"C: {counts['C']} ({percentages['C']:.2f}%)\t"
                f"W: {counts['W']} ({percentages['W']:.2f}%)\n"
            )
            # Write to output file
            out.write(output_line)

if __name__ == "__main__":
    input_fasta = "input_sequences.fasta"  # Replace with your input FASTA file
    output_file = "amino_acid_counts.tsv"
    
    process_fasta(input_fasta, output_file)
