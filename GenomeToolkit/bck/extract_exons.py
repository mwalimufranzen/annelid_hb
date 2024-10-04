from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re
import sys

# Function to extract sequence ID from the alignment file
def extract_sequence_id(alignment_file):
    with open(alignment_file, 'r') as file:
        content = file.read()
        match_id = re.search(r'Sequence ID:\s*(\S+)', content)
        if match_id:
            return match_id.group(1).strip()
        else:
            sys.stderr.write("No sequence ID found in the alignment file.\n")
            sys.exit(1)

# Function to download sequence from NCBI
def fetch_and_save_sequence(sequence_id, output_file):
    Entrez.email = 'your-email@example.com'  # Replace with your email
    try:
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

# Parse the alignment file to get ranges
def parse_alignment_file(alignment_file):
    ranges = []
    sequence_id = extract_sequence_id(alignment_file)

    with open(alignment_file, 'r') as file:
        content = file.read()

        # Extract ranges
        for match in re.finditer(r'Range \d+: (\d+) to (\d+)', content):
            start, end = map(int, match.groups())
            ranges.append((start, end))

    if not ranges:
        sys.stderr.write("No ranges found in the alignment file.\n")
        sys.exit(1)

    return sequence_id, ranges

# Extract exons and save to files
def extract_and_write_sequences(sequence_id, ranges, fasta_file):
    try:
        seq_record = SeqIO.read(fasta_file, "fasta")
        seq = str(seq_record.seq)
        sys.stdout.write(f"Sequence file '{fasta_file}' read successfully.\n")
    except Exception as e:
        sys.stderr.write(f"Error reading FASTA file: {e}\n")
        sys.exit(1)

    # Sort ranges based on their start positions for consistent processing
    ranges.sort(key=lambda x: min(x[0], x[1]))

    min_pos = float('inf')
    max_pos = float('-inf')
    records_written = []

    for i, (start, end) in enumerate(ranges):
        # Determine orientation based on the order of start and end
        if start > end:
            orientation = 'antisense'
            # For antisense, reverse complement the sequence from end to start
            # Adjust for 1-based to 0-based indexing
            range_seq = str(Seq(seq[end-1:start]).reverse_complement())
            # Preserve original start and end for naming
            original_start, original_end = end, start
        else:
            orientation = 'sense'
            range_seq = seq[start-1:end]
            original_start, original_end = start, end

        # Debug print statement for each range
        sys.stdout.write(f"Range {i+1}: start {original_start} end {original_end} orientation {orientation}\n")

        # Create exon filename
        exon_filename = f"exon{i+1}_{orientation}.fasta"
        exon_record = SeqRecord(Seq(range_seq),
                                id=f"{sequence_id}_{orientation}_range_{original_start}-{original_end}",
                                description=f"{orientation.capitalize()} Range {original_start}-{original_end}")
        with open(exon_filename, 'w') as file:
            SeqIO.write(exon_record, file, 'fasta')
            records_written.append(exon_filename)
        sys.stdout.write(f"Written {exon_filename}\n")

        # Create reverse complement file
        rc_filename = f"exon{i+1}rc_{orientation}.fasta"
        rc_seq = Seq(range_seq).reverse_complement()
        rc_record = SeqRecord(rc_seq,
                              id=f"{sequence_id}_rc_{orientation}_range_{original_start}-{original_end}",
                              description=f"Reverse complement {orientation.capitalize()} Range {original_start}-{original_end}")
        with open(rc_filename, 'w') as file:
            SeqIO.write(rc_record, file, 'fasta')
            records_written.append(rc_filename)
        sys.stdout.write(f"Written {rc_filename}\n")

        # Update min and max positions for flanking regions
        min_pos = min(min_pos, original_start)
        max_pos = max(max_pos, original_end)

    # Calculate flanking regions ensuring they don't exceed sequence boundaries
    flanking_start = max(min_pos - 350 - 1, 0)  # Subtract 1 for 0-based indexing
    flanking_end = min(max_pos + 350, len(seq))
    target_seq = seq[flanking_start:flanking_end]

    # Debug print statement for flanking regions
    sys.stdout.write(f"Flanking regions: start {flanking_start +1} end {flanking_end}\n")

    target_filename = "target.fasta"
    with open(target_filename, 'w') as file:
        target_record = SeqRecord(Seq(target_seq),
                                  id=f"{sequence_id}_truncated",
                                  description="Truncated sequence with additional flanking nucleotides")
        SeqIO.write(target_record, file, 'fasta')
        records_written.append(target_filename)
    sys.stdout.write(f"Written {target_filename}\n")

    # Inform the user about the files written
    sys.stdout.write("\nFiles written:\n")
    for filename in records_written:
        sys.stdout.write(f"{filename}\n")

# Main function
def main(alignment_file):
    sequence_id, ranges = parse_alignment_file(alignment_file)
    fasta_file = f"{sequence_id}.fasta"
    
    # Download the genomic sequence
    fetch_and_save_sequence(sequence_id, fasta_file)

    # Extract exons and write files
    extract_and_write_sequences(sequence_id, ranges, fasta_file)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Process alignment file and extract sequences.")
    parser.add_argument('alignment_file', help="File containing alignment ranges")
    args = parser.parse_args()

    main(args.alignment_file)
