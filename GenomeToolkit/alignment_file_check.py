from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re

def extract_sequence_id(alignment_file):
    """Extract sequence ID from the alignment file."""
    with open(alignment_file, 'r') as file:
        content = file.read()
        match_id = re.search(r'Sequence ID:\s*(\S+)', content)
        if match_id:
            return match_id.group(1).strip()
        else:
            raise ValueError("No sequence ID found in the alignment file.")

def parse_alignment_file(alignment_file):
    """Parse alignment file to get coding ranges and strand orientation."""
    ranges = []
    with open(alignment_file, 'r') as file:
        content = file.read()

        # Extract ranges and determine strand orientation
        for match in re.finditer(r'Range \d+: (\d+) to (\d+)', content):
            start, end = map(int, match.groups())
            # Extract only the first Sbjct line after each range
            sbjct_match = re.search(r'Sbjct\s+(\d+)\s+[^\n]+\s+(\d+)', content[match.end():])
            if sbjct_match:
                sbjct_start, sbjct_end = map(int, sbjct_match.groups())
                # Determine strand based on the first Sbjct number
                if sbjct_start == start:
                    orientation = 'sense'
                elif sbjct_start == end:
                    orientation = 'antisense'
                else:
                    raise ValueError(f"Unable to determine strand orientation for range {start}-{end}.")
                ranges.append((start, end, orientation))
            else:
                raise ValueError(f"No Sbjct line found for range {start}-{end}.")
    
    return ranges

def translate_sequence(seq, orientation):
    """Translate sequence based on orientation."""
    if orientation == 'antisense':
        # Reverse complement and translate
        return Seq(seq).reverse_complement().translate()
    else:
        # Translate sense strand
        return Seq(seq).translate()

def generate_open_reading_frames(seq, orientation):
    """Generate three open reading frames for the given orientation."""
    frames = []
    if orientation == 'antisense':
        seq = Seq(seq).reverse_complement()
    else:
        seq = Seq(seq)

    frames.append(seq.translate())
    frames.append(seq[1:].translate())
    frames.append(seq[2:].translate())
    
    return frames

def format_output(seq, max_line_length=60):
    """Format sequence output to a maximum line length."""
    return '\n'.join(seq[i:i+max_line_length] for i in range(0, len(seq), max_line_length))

def write_output(aa_sequences, orfs, output_file):
    """Write amino acid sequences and open reading frames to the output file."""
    with open(output_file, 'w') as out_file:
        # Write amino acid sequences for each range
        out_file.write("Translated Amino Acid Sequences for Ranges:\n")
        for i, aa_seq in enumerate(aa_sequences):
            out_file.write(f"\nRange {i+1}:\n")
            out_file.write(format_output(str(aa_seq)))
            out_file.write("\n")

        # Write open reading frames
        out_file.write("\nOpen Reading Frames:\n")
        for i, orf in enumerate(orfs):
            out_file.write(f"\nOpen reading frame {i+1}:\n")
            out_file.write(format_output(str(orf)))
            out_file.write("\n")

def main(alignment_file, cDNA_file, output_file):
    # Extract sequence ID from alignment file
    sequence_id = extract_sequence_id(alignment_file)
    fasta_file = f"{sequence_id}.fasta"

    # Parse alignment file to get ranges and strand orientation
    ranges = parse_alignment_file(alignment_file)

    # Read the cDNA or gene sequence from the fasta file
    seq_record = SeqIO.read(cDNA_file, "fasta")
    full_seq = str(seq_record.seq)

    # Extract and translate the coding regions based on the ranges
    aa_sequences = []
    orfs = []
    for start, end, orientation in ranges:
        coding_seq = full_seq[start-1:end]  # Extract coding sequence (1-based)
        aa_seq = translate_sequence(coding_seq, orientation)
        aa_sequences.append(aa_seq)

        # Generate open reading frames for the given sense
        orfs.extend(generate_open_reading_frames(coding_seq, orientation))

    # Write the amino acid sequences and open reading frames to the output file
    write_output(aa_sequences, orfs, output_file)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Translate tblastn alignment and compare to gene or cDNA.")
    parser.add_argument('alignment_file', help="File containing tblastn alignment ranges")
    parser.add_argument('cDNA_file', help="File containing the gene or cDNA sequence in FASTA format")
    parser.add_argument('output_file', help="Output file to write the results")
    args = parser.parse_args()

    main(args.alignment_file, args.cDNA_file, args.output_file)
