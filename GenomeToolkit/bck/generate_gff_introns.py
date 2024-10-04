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
    """Write the exon and intron positions to a GFF file."""
    with open(gff_file, "w") as gff:
        previous_exon_end = None

        for exon_name, (start, end) in exon_positions.items():
            if start is not None and end is not None:
                # Write exon information
                gff.write(f"{seqid}\tgenbank\t{exon_name}\t{start}\t{end}\t.\t+\t0\tattributes\n")
                
                # If there's a previous exon, write the intron between the exons
                if previous_exon_end is not None and previous_exon_end < start - 1:
                    intron_start = previous_exon_end + 1
                    intron_end = start - 1
                    intron_name = f"intron{int(exon_name[-1]) - 1}"  # Name the intron based on exons
                    gff.write(f"{seqid}\tgenbank\t{intron_name}\t{intron_start}\t{intron_end}\t.\t+\t0\tattributes\n")

                previous_exon_end = end

def main():
    # Read the target genomic sequence
    target_record = SeqIO.read("target.fasta", "fasta")
    target_seq = str(target_record.seq)
    seqid = target_record.id.split(".")[0]
    print(f"seqid from target: {seqid}")
    
    # Report the length of the target sequence
    print(f"Length of target sequence: {len(target_seq)}")

    # Read the exon sequences and report their lengths
    exon_files = ["exon1.fasta", "exon2.fasta", "exon3.fasta"]
    exon_positions = {}

    for i, exon_file in enumerate(exon_files, start=1):
        exon_record = SeqIO.read(exon_file, "fasta")
        exon_seq = str(exon_record.seq)
        
        print(f"Length of {exon_file}: {len(exon_seq)}")
        
        # Find start and end positions in the target sequence
        start, end = find_exon_positions(target_seq, exon_seq)
        
        if start is not None and end is not None:
            print(f"Exon {i} match found: start = {start}, end = {end}")
        else:
            print(f"Exon {i} match not found.")
        
        # Store the exon name and positions
        exon_positions[f"exon{i}"] = (start, end)

    # Write the results to a GFF file, including introns
    write_gff_file("tmp.gff", exon_positions, seqid)

if __name__ == "__main__":
    main()
