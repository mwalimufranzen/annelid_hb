from Bio import SeqIO

def find_exon_positions(target_seq, exon_seq):
    """Find the start and end positions of an exon within a target sequence."""
    start_seq = exon_seq[:8]  # First 8 nucleotides
    end_seq = exon_seq[-8:]   # Last 8 nucleotides
    
    print(f"Checking for start sequence: {start_seq}")
    start = target_seq.find(start_seq)
    
    if start != -1:
        # Now, check for the end sequence after the start
        print(f"Checking for end sequence: {end_seq}")
        end = target_seq.find(end_seq, start)
        if end != -1:
            # Adjust end position to include the full exon length
            end += len(end_seq)
            return start + 1, end  # 1-based offset
    return None, None

def write_gff_file(gff_file, exon_positions, seqid):
    """Write the exon positions to a GFF file."""
    with open(gff_file, "w") as gff:
        for exon_name, (start, end) in exon_positions.items():
            if start is not None and end is not None:
                gff.write(f"{seqid}\tgenbank\t{exon_name}\t{start}\t{end}\t.\t+\t0\tattributes\n")

def main():
    # Read the target genomic sequence
    target_record = SeqIO.read("target.fasta", "fasta")
    target_seq = str(target_record.seq)
    seqid = target_record.id.split(".")[0]
    print(f"seqid from target: {seqid}")
    
    # Report the length of the target sequence
    print(f"Length of target sequence: {len(target_seq)}")

    # Initialize dictionary to store exon positions
    exon_positions = {}

    # Process exon1
    exon1_record = SeqIO.read("exon1.fasta", "fasta")
    exon1_seq = str(exon1_record.seq)
    print(f"Length of exon1: {len(exon1_seq)}")
    start1, end1 = find_exon_positions(target_seq, exon1_seq)
    if start1 is not None and end1 is not None:
        print(f"Exon1 match found: start = {start1}, end = {end1}")
    else:
        print(f"Exon1 match not found.")
    exon_positions["exon1"] = (start1, end1)

    # Process exon2
    exon2_record = SeqIO.read("exon2.fasta", "fasta")
    exon2_seq = str(exon2_record.seq)
    print(f"Length of exon2: {len(exon2_seq)}")
    start2, end2 = find_exon_positions(target_seq, exon2_seq)
    if start2 is not None and end2 is not None:
        print(f"Exon2 match found: start = {start2}, end = {end2}")
    else:
        print(f"Exon2 match not found.")
    exon_positions["exon2"] = (start2, end2)

    # Write the results to a GFF file
    write_gff_file("tmp.gff", exon_positions, seqid)

if __name__ == "__main__":
    main()
