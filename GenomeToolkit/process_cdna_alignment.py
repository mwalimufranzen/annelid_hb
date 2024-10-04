from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re
import sys
import requests
import warnings
from Bio import BiopythonWarning

# Ignore the Biopython warnings about partial codons for now
warnings.simplefilter('ignore', BiopythonWarning)

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

# Function to download the sequence file from NCBI based on the sequence ID
def download_sequence(sequence_id, output_file):
    url = f"https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id={sequence_id}&db=nuccore&report=fasta"
    try:
        response = requests.get(url)
        if response.status_code == 200:
            with open(output_file, 'w') as file:
                file.write(response.text)
            sys.stdout.write(f"Sequence {sequence_id} downloaded successfully.\n")
        else:
            sys.stderr.write(f"Failed to download sequence for ID {sequence_id}. HTTP Status: {response.status_code}\n")
            sys.exit(1)
    except Exception as e:
        sys.stderr.write(f"Error downloading sequence: {e}\n")
        sys.exit(1)

# Parse the alignment file to get ranges and strand orientation
def parse_alignment_file(alignment_file):
    ranges = []
    sequence_id = extract_sequence_id(alignment_file)

    with open(alignment_file, 'r') as file:
        content = file.read()

        # Extract ranges and determine strand orientation based on the first occurrence of Sbjct line
        for match in re.finditer(r'Range \d+: (\d+) to (\d+)', content):
            start, end = map(int, match.groups())
            # Extract only the first Sbjct line after each range using a flexible regex
            sbjct_match = re.search(r'Sbjct\s+(\d+)\s+[^\n]+\s+(\d+)', content[match.end():])
            if sbjct_match:
                sbjct_start, sbjct_end = map(int, sbjct_match.groups())
                # Determine strand based on the first Sbjct number
                if sbjct_start == start:
                    orientation = 'sense'
                elif sbjct_start == end:
                    orientation = 'antisense'
                else:
                    sys.stderr.write(f"Unable to determine strand orientation for range {start}-{end}.\n")
                    sys.stderr.write(f"Sbjct start: {sbjct_start}, Sbjct end: {sbjct_end}, expected start: {start}, expected end: {end}\n")
                    sys.exit(1)
                ranges.append((start, end, orientation))
                sys.stdout.write(f"Range {len(ranges)}: start {start}, end {end}, {orientation}\n")
            else:
                sys.stderr.write(f"No Sbjct line found for range {start}-{end}. Trying to proceed.\n")
                ranges.append((start, end, 'unknown'))  # Allow it to proceed even without orientation

    # Sort ranges by start position to handle out-of-order ranges
    ranges = sorted(ranges, key=lambda x: x[0])

    if not ranges:
        sys.stderr.write("No ranges found in the alignment file.\n")
        sys.exit(1)

    return sequence_id, ranges

# Extract exons, introns, and write sequences to output files
def extract_exons_introns(alignment_file, ranges, fasta_file):
    try:
        seq_record = SeqIO.read(fasta_file, "fasta")
        seq = str(seq_record.seq)
        sys.stdout.write(f"Sequence file '{fasta_file}' read successfully.\n")
    except Exception as e:
        sys.stderr.write(f"Error reading FASTA file: {e}\n")
        sys.exit(1)

    min_pos = float('inf')
    max_pos = float('-inf')

    exon_positions = []
    with open(alignment_file, 'r') as aln_file:
        alignment_data = aln_file.read()

    # Iterate over ranges to extract exons and determine orientation
    for i, (start, end, orientation) in enumerate(ranges):
        exon_positions.append((start, end))  # Save positions for intron extraction

        if orientation == 'antisense':
            # Antisense: extract from start to end and reverse complement properly
            range_seq = str(Seq(seq[start-1:end]).reverse_complement())
            exon_filename = f"exon{i+1}_antisense.fasta"
        elif orientation == 'sense':
            # Sense strand
            range_seq = seq[start-1:end]
            exon_filename = f"exon{i+1}_sense.fasta"
        else:
            sys.stderr.write(f"Unknown strand orientation for range {start}-{end}. Skipping.\n")
            continue

        # Write reverse complement sequence to exon files
        with open(exon_filename, 'w') as exon_file:
            SeqIO.write(SeqRecord(Seq(range_seq), id=f"exon{i+1}", description=f"Range {start}-{end}"), exon_file, 'fasta')

        sys.stdout.write(f"Written exon file for range {start}-{end}.\n")

        # Track the min and max positions for target flanking
        min_pos = min(min_pos, start)
        max_pos = max(max_pos, end)

    # Write target file with 350 bp flanking regions
    flanking_start = max(min_pos - 350, 0)
    flanking_end = min(max_pos + 350, len(seq))
    target_seq = seq[flanking_start:flanking_end]

    target_filename = "target.fasta"
    with open(target_filename, 'w') as file:
        target_record = SeqRecord(Seq(target_seq),
                                  id="target",
                                  description=f"Target sequence with 350bp flanking regions")
        SeqIO.write(target_record, file, 'fasta')

    sys.stdout.write(f"Target sequence written to {target_filename}.\n")

    # Extract and write introns based on exon positions (always in sense orientation)
    if len(exon_positions) > 1:
        for j in range(1, len(exon_positions)):
            intron_start = exon_positions[j-1][1] + 1
            intron_end = exon_positions[j][0] - 1
            if intron_start < intron_end:
                intron_seq = seq[intron_start-1:intron_end]
                intron_filename = f"intron{j}.fasta"
                with open(intron_filename, 'w') as intron_file:
                    SeqIO.write(SeqRecord(Seq(intron_seq), id=f"intron{j}", description=f"Intron {j} between exons"), intron_file, 'fasta')
                sys.stdout.write(f"Intron {j} written to {intron_filename}.\n")

# Main function
def main(alignment_file):
    sequence_id, ranges = parse_alignment_file(alignment_file)
    fasta_file = f"{sequence_id}.fasta"

    # Download the FASTA file if it doesn't exist
    download_sequence(sequence_id, fasta_file)

    # Extract exons, introns, and write sequences to files
    extract_exons_introns(alignment_file, ranges, fasta_file)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Process alignment file and extract nucleotide sequences.")
    parser.add_argument('alignment_file', help="File containing alignment ranges")
    args = parser.parse_args()

    main(args.alignment_file)
