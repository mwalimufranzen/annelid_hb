import re
from Bio import SeqIO

FLANKING_LENGTH = 350  # Flanking sequence in target.fasta
SEARCH_RANGE = 10  # Number of nucleotides to search around the boundary for the correct splice site

def parse_alignment_file(alignment_file):
    """Parse the alignment.txt file to extract the Protein ID and determine sense or antisense."""
    with open(alignment_file, 'r') as f:
        content = f.read()

    # Extract the Protein ID
    protein_id_match = re.search(r'Query #1: tr\|(\w+)\|', content)
    if protein_id_match:
        protein_id = protein_id_match.group(1)
        print(f"Protein ID found: {protein_id}")
    else:
        raise ValueError("Protein ID not found in alignment.txt")

    # Extract the Range and Sbjct values to determine sense or antisense
    range_match = re.search(r'Range \d+: (\d+) to (\d+)', content)
    sbjct_match = re.search(r'Sbjct\s+(\d+)\s+[^\n]+\s+(\d+)', content)
    
    if range_match and sbjct_match:
        start, end = map(int, range_match.groups())
        sbjct_start, sbjct_end = map(int, sbjct_match.groups())
        
        if sbjct_start == start:
            orientation = "+"
        elif sbjct_start == end:
            orientation = "-"
        else:
            raise ValueError("Unable to determine sense or antisense orientation.")
        
        print(f"Orientation determined: {'sense' if orientation == '+' else 'antisense'}")
        return protein_id, orientation
    else:
        raise ValueError("Range or Sbjct information not found in alignment.txt")

def read_gff_file(gff_file):
    """Read the original GFF file to extract exon and intron positions."""
    exon_positions = []
    sequence_id = None

    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if fields[2].startswith("exon"):
                exon_name = fields[2]
                start = int(fields[3])
                end = int(fields[4])
                orientation = fields[6]
                sequence_id = fields[0]
                exon_positions.append((exon_name, start, end, orientation))
    
    return sequence_id, exon_positions

def adjust_intron_boundaries(seq, start, end, orientation):
    """Adjust the intron boundaries based on the closest splice site (GT-AG or CT-AC for antisense).
       The search covers 6 nucleotides on either side of the original boundary."""
    SEARCH_SPAN = 6  # We will look 6 nucleotides on either side of the boundary

    # Initialize variables to track the closest positions
    closest_start = start
    closest_end = end
    closest_start_distance = SEARCH_SPAN + 1  # Set larger than max search span
    closest_end_distance = SEARCH_SPAN + 1  # Set larger than max search span

    if orientation == '+':
        # Search for GT-AG splice sites around the boundary
        for i in range(-SEARCH_SPAN, SEARCH_SPAN + 1):
            # Check GT for the start
            if seq[start + i:start + i + 2] == 'GT':
                if abs(i) < closest_start_distance:
                    closest_start = start + i
                    closest_start_distance = abs(i)

            # Check AG for the end
            if seq[end + i:end + i + 2] == 'AG':
                if abs(i) < closest_end_distance:
                    closest_end = end + i
                    closest_end_distance = abs(i)
    else:
        # For antisense, look for CT-AC splice sites
        for i in range(-SEARCH_SPAN, SEARCH_SPAN + 1):
            # Check CT for the start
            if seq[start + i:start + i + 2] == 'CT':
                if abs(i) < closest_start_distance:
                    closest_start = start + i
                    closest_start_distance = abs(i)

            # Check AC for the end
            if seq[end + i:end + i + 2] == 'AC':
                if abs(i) < closest_end_distance:
                    closest_end = end + i
                    closest_end_distance = abs(i)

    # Adjust positions to 1-based indexing for GFF output
    return closest_start + 1, closest_end + 2

def adjust_exon_intron_boundaries(exon_positions, seq):
    """Adjust the exon and intron boundaries while ensuring exon 3 start and exon 1 end are fixed."""
    intron_positions = []

    # Iterate through exons to calculate introns between them
    for idx, (exon_name, start, end, orientation) in enumerate(exon_positions):
        if idx > 0:
            if idx == 2:  # Adjust Exon 3 and Intron 2
                # Intron 2: starts after exon 3 and ends before exon 2
                intron_start = exon_positions[2][2] + 1  # End of exon 3 + 1
                intron_end = exon_positions[1][1] - 1  # Start of exon 2 - 1
                intron_start, intron_end = adjust_intron_boundaries(seq, intron_start, intron_end, orientation)

                # Fix Exon 3 start and adjust Exon 2 start
                exon_positions[2] = (exon_positions[2][0], exon_positions[2][1], intron_start - 1, orientation)  # Exon 3 end
                exon_positions[1] = (exon_positions[1][0], intron_end + 1, exon_positions[1][2], orientation)  # Exon 2 start

                intron_positions.append((f"intron2", intron_start, intron_end, orientation))

            if idx == 1:  # Adjust Exon 1 and Intron 1
                # Intron 1: starts after exon 2 and ends before exon 1
                intron_start = exon_positions[1][2] + 1  # End of exon 2 + 1
                intron_end = exon_positions[0][1] - 1  # Start of exon 1 - 1
                intron_start, intron_end = adjust_intron_boundaries(seq, intron_start, intron_end, orientation)

                # Fix Exon 1 start and adjust Exon 2 end
                exon_positions[1] = (exon_positions[1][0], exon_positions[1][1], intron_start - 1, orientation)  # Exon 2 end
                exon_positions[0] = (exon_positions[0][0], intron_end + 1, exon_positions[0][2], orientation)  # Exon 1 start

                intron_positions.append((f"intron1", intron_start, intron_end, orientation))

            # Debugging output
            print(f"Adjusted intron {idx}: start {intron_start}, end {intron_end}")
            print(f"Adjusted exon {idx+1}: start {exon_positions[idx][1]}, end {exon_positions[idx][2]}")

    return intron_positions


def check_start_stop_codons(exon_positions, seq, orientation):
    """Check for start and stop codons in the exons."""
    if orientation == '+':
        start_codon = seq[exon_positions[0][1] - 1:exon_positions[0][1] + 2]
        if start_codon == 'ATG':
            print("Start codon (ATG) found at the start of exon 1.")
        
        stop_codon = seq[exon_positions[-1][2]:exon_positions[-1][2] + 3]
        if stop_codon in ['TAA', 'TAG', 'TGA']:
            print(f"Stop codon ({stop_codon}) found at the end of exon 3.")
    else:
        start_codon = seq[exon_positions[-1][2] - 3:exon_positions[-1][2]]
        if start_codon == 'TAC':
            print("Start codon (TAC) found at the end of exon 3.")

        stop_codon = seq[exon_positions[0][1] - 4:exon_positions[0][1] - 1]
        if stop_codon in ['TTA', 'CTA', 'TCA']:
            print(f"Stop codon ({stop_codon}) found before exon 1.")

def write_gff_file(gff_file, exon_positions, intron_positions, sequence_id, protein_id):
    """Write the exon and intron positions to a new GFF file."""
    with open(gff_file, "w") as gff:
        for exon_name, start, end, orientation in exon_positions:
            gff.write(f"{sequence_id}\tgenbank\t{exon_name}\t{start}\t{end}\t.\t{orientation}\t.\tProtein ID: {protein_id}; exon_number={exon_name[-1]}\n")
        
        for intron_name, start, end, orientation in intron_positions:
            gff.write(f"{sequence_id}\tgenbank\t{intron_name}\t{start}\t{end}\t.\t.\t.\tProtein ID: {protein_id}; intron_number={intron_name[-1]}\n")


def excise_exons_to_fasta(exon_positions, target_seq, protein_id):
    """Excise the exons from the genomic data and write them to exon#.fasta files using the correct adjusted limits."""
    for idx, (exon_name, start, end, orientation) in enumerate(exon_positions):
        # Use the adjusted exon positions (from the GFF file logic)
        exon_seq = target_seq[start - 1:end]  # Adjust for 0-based indexing

        # Write the sense file (always the same for both orientations)
        sense_file_name = f"{exon_name}_sense.fasta"
        with open(sense_file_name, "w") as exon_file:
            exon_file.write(f">{exon_name}_sense {start}-{end}\n{exon_seq}\n")
        print(f"Exon {idx + 1} written to {sense_file_name}: start={start}, end={end}")

        # If the strand is antisense, also write the reverse complement
        if orientation == '-':
            antisense_file_name = f"{exon_name}_antisense.fasta"
            reverse_complement_seq = Seq(exon_seq).reverse_complement()  # Convert to Seq and get reverse complement
            with open(antisense_file_name, "w") as exon_file:
                exon_file.write(f">{exon_name}_antisense {start}-{end}\n{reverse_complement_seq}\n")
            print(f"Reverse complement exon written to {antisense_file_name}: start={start}, end={end}")

def write_intron_files(exon_positions, target_seq):
    """Extract and write introns based on exon positions with sense/antisense designations using adjusted limits."""
    if len(exon_positions) > 1:
        for j in range(1, len(exon_positions)):
            # Calculate intron positions using the correct boundaries
            intron_start = exon_positions[j-1][2] + 1  # End of previous exon +1
            intron_end = exon_positions[j][1] - 1  # Start of next exon -1

            print(f"Intron {j} range: start {intron_start}, end {intron_end}")

            if intron_start < intron_end:
                intron_seq = target_seq[intron_start-1:intron_end]  # Adjust for 0-based indexing
                orientation = exon_positions[j-1][3]  # Use the orientation of the previous exon

                # Write the sense file (same for both orientations)
                sense_intron_filename = f"intron{j}_sense.fasta"
                with open(sense_intron_filename, 'w') as intron_file:
                    intron_file.write(f">intron{j}_sense {intron_start}-{intron_end}\n{intron_seq}\n")
                print(f"Intron {j} written to {sense_intron_filename}: start={intron_start}, end={intron_end}")

                # If the orientation is antisense, also write the reverse complement
                if orientation == '-':
                    antisense_intron_filename = f"intron{j}_antisense.fasta"
                    reverse_complement_intron = Seq(intron_seq).reverse_complement()  # Convert to Seq and get reverse complement
                    with open(antisense_intron_filename, 'w') as intron_file:
                        intron_file.write(f">intron{j}_antisense {intron_start}-{intron_end}\n{reverse_complement_intron}\n")
                    print(f"Reverse complement intron written to {antisense_intron_filename}: start={intron_start}, end={intron_end}")


def main():
    alignment_file = "alignment.txt"
    target_file = "target.fasta"

    # Step 1: Retrieve Protein ID and orientation (sense or antisense) from alignment.txt
    protein_id, orientation = parse_alignment_file(alignment_file)
    
    # Step 2: Read the original GFF file for this Protein ID
    original_gff_file = f"{protein_id}.gff"
    sequence_id, exon_positions = read_gff_file(original_gff_file)

    # Step 3: Read the target sequence from target.fasta
    target_record = SeqIO.read(target_file, "fasta")
    full_seq = str(target_record.seq)

    # Step 4: Adjust the intron-exon boundaries and ensure they are sequential
    intron_positions = adjust_exon_intron_boundaries(exon_positions, full_seq)

    # Step 5: Check for start and stop codons
    check_start_stop_codons(exon_positions, full_seq, orientation)

    # Step 6: Write new GFF file
    output_gff_file = f"{protein_id}_chk.gff"
    write_gff_file(output_gff_file, exon_positions, intron_positions, sequence_id, protein_id)

    # Step 7: Write exon and intron FASTA files
    excise_exons_introns_to_fasta(exon_positions, intron_positions, full_seq)

    print(f"Adjusted GFF file '{output_gff_file}' and exon/intron FASTA files have been created successfully.")

if __name__ == "__main__":
    main()
