from Bio import SeqIO

def count_amino_acids(seq):
    """Count the specific amino acids in a sequence."""
    counts = {
        'M': seq.count('M'),
        'Y': seq.count('Y'),
        'H': seq.count('H'),
        'C': seq.count('C'),
        'W': seq.count('W')
    }
    return counts

def process_fasta(input_fasta, output_file):
    """Process the FASTA file and write the amino acid counts to the output file."""
    with open(output_file, 'w') as out:
        for record in SeqIO.parse(input_fasta, "fasta"):
            # Extract the accession number after 'tr|'
            accession = record.id.split('|')[1] if '|' in record.id else record.id
            # Count the amino acids
            counts = count_amino_acids(str(record.seq))
            # Format the output line
            output_line = f"{accession}\tM: {counts['M']}\tY: {counts['Y']}\tH: {counts['H']}\tC: {counts['C']}\tW: {counts['W']}\n"
            # Write to output file
            out.write(output_line)

if __name__ == "__main__":
    input_fasta = "input_sequences.fasta"  # Replace with your input FASTA file
    output_file = "amino_acid_counts.tsv"
    
    process_fasta(input_fasta, output_file)
