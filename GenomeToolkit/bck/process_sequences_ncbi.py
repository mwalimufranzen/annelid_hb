from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re
import sys
import os

def fetch_sequence_from_ncbi(sequence_id):
    Entrez.email = 'franzen@ncsu.edu'  # Replace with your email
    try:
        sys.stderr.write(f"Fetching sequence for ID: {sequence_id}\n")
        handle = Entrez.efetch(db="nuccore", id=sequence_id, rettype="fasta", retmode="text")
        return handle.read()
    except Exception as e:
        sys.stderr.write(f"Error fetching sequence from NCBI: {e}\n")
        sys.exit(1)

def extract_ranges_and_id(alignment_text):
    ranges = []
    sequence_id = None
    organism = None
    lines = alignment_text.strip().split('\n')
    
    for line in lines:
        if line.startswith('Sequence ID:'):
            sequence_id = line.split(':')[1].strip()
        elif line.startswith('Range'):
            match = re.match(r'Range \d+: (\d+) to (\d+)', line)
            if match:
                start, end = map(int, match.groups())
                ranges.append((start, end))
        elif 'OS=' in line:
            organism = line.split('OS=')[1].strip()
    
    if not organism:
        sys.stderr.write("No organism information found in the alignment file.\n")
    
    return sequence_id, organism, ranges

def sanitize_filename(name, max_length=100):
    # Create a sanitized filename
    sanitized = re.sub(r'[^\w\-]', '_', name)  # Replace non-alphanumeric characters
    if len(sanitized) > max_length:
        sanitized = sanitized[:max_length]  # Truncate if too long
    return sanitized

def extract_and_write_sequences(sequence_id, ranges, seq, output_fasta_file, truncated_file):
    min_pos = float('inf')
    max_pos = float('-inf')
    
    records = []
    
    # Use shortened ID for filenames
    short_id = sanitize_filename(sequence_id)
    
    sys.stdout.write(f"Sequence ID length: {len(sequence_id)}\n")
    sys.stdout.write(f"Ranges found in the alignment file:\n")
    for start, end in ranges:
        sys.stdout.write(f"Range: {start} to {end}\n")
    
    for i, (start, end) in enumerate(ranges):
        if start > end:
            sys.stderr.write(f"Warning: Start position {start} is greater than end position {end} in range {i+1}. Skipping this range.\n")
            continue
        
        # Ensure that the start and end are within bounds
        start = max(start - 1, 0)
        end = min(end, len(seq))
        
        range_seq = seq[start:end]
        
        min_pos = min(min_pos, start + 1)  # Convert to 1-based indexing for reporting
        max_pos = max(max_pos, end)
        
        record = SeqRecord(Seq(range_seq),
                           id=f"{short_id}_range_{start+1}-{end}",
                           description=f"Range {start+1}-{end}")
        records.append(record)
    
    try:
        with open(output_fasta_file, 'w') as fasta_out:
            SeqIO.write(records, fasta_out, 'fasta')
        sys.stderr.write(f"Multifasta file written to: {output_fasta_file}\n")
    except Exception as e:
        sys.stderr.write(f"Error writing multifasta file: {e}\n")
        sys.exit(1)
    
    truncated_start = max(min_pos - 350, 1)
    truncated_end = min(max_pos + 350, len(seq))
    
    truncated_seq = seq[truncated_start-1:truncated_end]
    
    try:
        with open(truncated_file, 'w') as fasta_out:
            truncated_record = SeqRecord(Seq(truncated_seq),
                                         id=f"{short_id}_truncated",
                                         description="Truncated sequence")
            SeqIO.write(truncated_record, fasta_out, 'fasta')
        sys.stderr.write(f"Truncated file written to: {truncated_file}\n")
    except Exception as e:
        sys.stderr.write(f"Error writing truncated file: {e}\n")
        sys.exit(1)

def main(alignment_file, output_fasta_file, truncated_file):
    try:
        with open(alignment_file, 'r') as align_file:
            alignment_text = align_file.read()
    except FileNotFoundError:
        sys.stderr.write(f"Alignment file '{alignment_file}' not found.\n")
        sys.exit(1)
    
    sequence_id, organism, ranges = extract_ranges_and_id(alignment_text)
    
    if not sequence_id:
        sys.stderr.write("No valid sequence ID found in the alignment file.\n")
        sys.exit(1)
    
    if not organism:
        sys.stderr.write("No organism information found in the alignment file.\n")
    
    fasta_data = fetch_sequence_from_ncbi(sequence_id)
    
    try:
        # Read the sequence data from FASTA format
        seq_record = SeqIO.read(fasta_data, "fasta")
        seq = str(seq_record.seq)
    except Exception as e:
        sys.stderr.write(f"Error parsing FASTA data: {e}\n")
        sys.exit(1)
    
    extract_and_write_sequences(sequence_id, ranges, seq, output_fasta_file, truncated_file)

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Process sequence alignments and extract ranges.")
    parser.add_argument('alignment_file', help="File containing alignment ranges")
    parser.add_argument('output_fasta_file', help="Output file for the multifasta sequences")
    parser.add_argument('truncated_file', help="Output file for the truncated sequence")
    
    args = parser.parse_args()
    
    main(args.alignment_file, args.output_fasta_file, args.truncated_file)
