from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re
import sys

def parse_alignment_file(alignment_file):
    ranges = []
    sequence_id = None
    with open(alignment_file, 'r') as file:
        content = file.read()
        # Extract sequence ID
        match_id = re.search(r'Sequence ID:\s*(\S+)', content)
        if match_id:
            sequence_id = match_id.group(1).strip()
        else:
            sys.stderr.write("No sequence ID found in the alignment file.\n")
            sys.exit(1)
        
        # Extract ranges
        for match in re.finditer(r'Range \d+: (\d+) to (\d+)', content):
            start, end = map(int, match.groups())
            ranges.append((start, end))
    
    if not ranges:
        sys.stderr.write("No ranges found in the alignment file.\n")
        sys.exit(1)

    return sequence_id, ranges

def extract_and_write_sequences(sequence_id, ranges, fasta_file):
    try:
        # Read the genomic sequence
        seq_record = SeqIO.read(fasta_file, "fasta")
        seq = str(seq_record.seq)
    except Exception as e:
        sys.stderr.write(f"Error reading FASTA file: {e}\n")
        sys.exit(1)

    # Sort ranges by start position
    ranges.sort()

    min_pos = float('inf')
    max_pos = float('-inf')
    records = []
    
    for i, (start, end) in enumerate(ranges):
        if start > end:
            sys.stderr.write(f"Warning: Start position {start} is greater than end position {end} in range {i+1}. Skipping this range.\n")
            continue
        
        # Adjust indices for 1-based indexing and extract sequence
        start -= 1
        end -= 1
        range_seq = seq[start:end + 1]
        
        # Reverse the sequence for antisense (reverse complement)
        antisense_seq = str(Seq(range_seq).reverse_complement())
        
        # Save antisense (reverse complement) version
        antisense_record = SeqRecord(Seq(antisense_seq),
                           id=f"{sequence_id}_antisense_range_{start+1}-{end+1}",
                           description=f"Antisense Range {start+1}-{end+1}")
        with open(f"exon{i+1}_antisense.fasta", 'w') as file:
            SeqIO.write(antisense_record, file, 'fasta')
        
        # Save sense (transcribed) version
        sense_seq = str(Seq(range_seq))
        sense_record = SeqRecord(Seq(sense_seq),
                           id=f"{sequence_id}_sense_range_{start+1}-{end+1}",
                           description=f"Sense Range {start+1}-{end+1}")
        with open(f"exon{i+1}_sense.fasta", 'w') as file:
            SeqIO.write(sense_record, file, 'fasta')
        
        min_pos = min(min_pos, start + 1)  # Convert to 1-based indexing
        max_pos = max(max_pos, end + 1)
    
    # Write target file with additional flanking nucleotides
    truncated_start = max(min_pos - 350, 0)
    truncated_end = min(max_pos + 350, len(seq))
    target_seq = seq[truncated_start:truncated_end]
    
    with open("target.fasta", 'w') as file:
        target_record = SeqRecord(Seq(target_seq),
                                  id=f"{sequence_id}_truncated",
                                  description="Truncated sequence with additional flanking nucleotides")
        SeqIO.write(target_record, file, 'fasta')

def main(alignment_file, fasta_file):
    sequence_id, ranges = parse_alignment_file(alignment_file)
    
    # Print ranges to monitor
    sys.stdout.write("Ranges found in the alignment file:\n")
    for start, end in ranges:
        sys.stdout.write(f"Range: {start} to {end}\n")
    
    extract_and_write_sequences(sequence_id, ranges, fasta_file)

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Process alignment file and extract sequences.")
    parser.add_argument('alignment_file', help="File containing alignment ranges")
    parser.add_argument('fasta_file', help="FASTA file with the genomic sequence")
    
    args = parser.parse_args()
    
    main(args.alignment_file, args.fasta_file)
