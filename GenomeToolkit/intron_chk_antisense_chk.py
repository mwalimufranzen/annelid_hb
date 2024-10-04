import re
from Bio import SeqIO
from Bio.Seq import Seq  # For reverse complement

FLANKING_LENGTH = 350  # Flanking sequence in target.fasta
SEARCH_RANGE = 20  # Number of nucleotides to search around the boundary for the correct splice site

def parse_alignment_file(alignment_file):
    """Parse the alignment.txt file to extract the Protein ID and determine sense or antisense."""
    with open(alignment_file, 'r') as f:
        content = f.read()

    # Extract the Protein ID
    #protein_id_match = re.search(r'Query #1: tr\|(\w+)\|', content)
    protein_id_match = re.search(r'Protein ID: (\w+)', content)
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
    
    # Report the number of exons read
    print(f"Number of exon positions read from the GFF file: {len(exon_positions)}")
    
    return sequence_id, exon_positions


def adjust_intron_boundaries(exon_positions, seq):
    """Adjust intron boundaries based on GT-AG (sense) or CT-AC (antisense) splice sites without changing exons."""
    intron_positions = []
    is_antisense = exon_positions[0][3] == '-'

    if is_antisense:
        print(f"Orientation {exon_positions[0][3]}")
        # Iterate in reverse order for antisense
        for idx in range(1, len(exon_positions)):
            # Reverse indexing for antisense
            exon_idx = len(exon_positions) - idx
            prev_exon = exon_positions[exon_idx - 1]
            current_exon = exon_positions[exon_idx]

            # Define initial intron boundaries according to input exon positions
            intron_start = current_exon[2] + 1  # End of the current exon + 1
            intron_end = prev_exon[1] - 1  # Start of the previous exon - 1

            # Adjust the intron boundaries based on the sequence
            intron_start, intron_end = find_splice_sites(seq, intron_start, intron_end, current_exon[3])

            # Append the adjusted intron to the list of intron positions
            intron_positions.append((f"intron{len(exon_positions) - idx}", intron_start, intron_end, current_exon[3]))

            # Debugging output
            print(f"Intron {len(exon_positions) - idx}: Adjusted Start={intron_start}, Adjusted End={intron_end}")

    else:  
        # Iterate through exons to calculate introns between them in the sense direction
        #NOT DEBUGGED FOR SENSE 
        for idx in range(1, len(exon_positions)):
            prev_exon = exon_positions[idx - 1]
            current_exon = exon_positions[idx]

            intron_start = prev_exon[2] + 1  # End of the previous exon + 1
            intron_end = current_exon[1] - 1  # Start of the current exon - 1

            # Adjust the intron boundaries based on the sequence
            intron_start, intron_end = find_splice_sites(seq, intron_start, intron_end, current_exon[3])

            # Append the adjusted intron to the list of intron positions
            intron_positions.append((f"intron{idx}", intron_start, intron_end, current_exon[3]))

            # Debugging output
            print(f"Intron {idx}: Adjusted Start={intron_start}, Adjusted End={intron_end}")

    return intron_positions

def adjust_exon_boundaries(exon_positions, intron_positions, is_antisense):
    """Adjust exon boundaries based on the new intron boundaries."""

    if is_antisense:
        # Antisense: Exon 4 (start4, end4) > intron 3 (end4 + 1, start3 -1) > Exon 3 (start3, end3) ...
        for idx in range(len(exon_positions)):
            exon_idx = len(exon_positions) - idx -1  # Reverse index for antisense
            
            # intron indices are the reverse of exon
            if(idx < len(intron_positions)):
                current_intron_start, current_intron_end = intron_positions[idx][1], intron_positions[idx][2]
                #print(f"Index {exon_idx + 1}: Current intron start ={current_intron_start}, Current intron end={current_intron_end}")
            if(idx > 0):
                previous_intron_start, previous_intron_end = intron_positions[idx - 1][1], intron_positions[idx - 1][2]
                #print(f"Index {exon_idx + 1}: Previous intron start ={previous_intron_start}, Previous intron end={previous_intron_end}")
            
            # Handle exon 4 (highest exon)
            if exon_idx == len(exon_positions) - 1:  # First exon in reverse (highest exon)
                exon_positions[exon_idx] = (
                    exon_positions[exon_idx][0],
                    exon_positions[exon_idx][1],  # Start stays the same for the highest exon
                    current_intron_start - 1,  # Adjust the end of the current exon based on previous intron
                    exon_positions[exon_idx][3]
                )
                print(f"Exon {exon_idx + 1}: Adjusted Start={exon_positions[exon_idx][1]}, Adjusted End={current_intron_start - 1}")
            # Handle internal exons
            elif exon_idx > 0:
                exon_positions[exon_idx] = (
                    exon_positions[exon_idx][0],
                    previous_intron_end + 1,  # Adjust the start of the current exon based on current intron
                    current_intron_start - 1,  # Adjust the end of the current exon based on current intron
                    exon_positions[exon_idx][3]
                )
                print(f"Exon {exon_idx + 1}: Adjusted Start={previous_intron_end + 1}, Adjusted End={current_intron_start - 1}")
            # Handle exon 1 (lowest exon in reverse)
            elif exon_idx == 0:
                exon_positions[exon_idx] = (
                    exon_positions[exon_idx][0],
                    previous_intron_end + 1,  # Adjust the start of the current exon based on current intron
                    exon_positions[exon_idx][2],  # End stays the same for the lowest exon
                    exon_positions[exon_idx][3]
                )
                print(f"Exon {exon_idx + 1}: Adjusted Start={previous_intron_end + 1}, Adjusted End={exon_positions[exon_idx][2]}")

    else:
        # Sense: Exon 4 (start4, end4) > intron 3 (start4 -1, end3 + 1) > Exon 3 (start3, end3) ...
        for idx in range(len(exon_positions)):
            if(idx < len(intron_positions)):
                current_intron_start, current_intron_end = intron_positions[idx][1], intron_positions[idx][2]
                print(f"Index {exon_idx + 1}: Current intron start ={current_intron_start}, Current intron end={current_intron_end}")
            if(idx > 0):
                previous_intron_start, previous_intron_end = intron_positions[idx - 1][1], intron_positions[idx - 1][2]
                print(f"Index {exon_idx + 1}: Previous intron start ={previous_intron_start}, Previous intron end={previous_intron_end}") 

            # Handle exon 1 (first exon in forward)
            if idx == 0:  # First exon
                
                exon_positions[idx] = (
                    exon_positions[idx][0],
                    exon_positions[idx][1],  # Start stays the same for the first exon
                    current_intron_start - 1,  # Adjust the end of the current exon based on current intron
                    exon_positions[idx][3]
                )
                print(f"Exon {idx + 1}: Adjusted Start={exon_positions[idx][1]}, Adjusted End={current_intron_start - 1}")
            # Handle internal exons
            elif idx < len(exon_positions):
                exon_positions[idx] = (
                    exon_positions[idx][0],
                    previous_intron_end + 1,  # Adjust the start of the current exon based on current intron
                    current_intron_start - 1,  # Adjust the end of the current exon based on current intron
                    exon_positions[idx][3]
                )
                print(f"Exon {idx + 1}: Adjusted Start={current_intron_end + 1}, Adjusted End={current_intron_start - 1}")
            elif idx == len(exon_positions):
                exon_positions[idx] = (
                    exon_positions[idx][0],
                    previous_intron_end + 1,  # Adjust the start of the current exon based on current intron
                    exon_positions[idx][1],  # Adjust the end of the current exon based on current intron
                    exon_positions[idx][3]
                )
                print(f"Exon {idx + 1}: Adjusted Start={current_intron_end + 1}, Adjusted End={exon_positions[idx][1]}")
                # Adjust the previous exon start

    return exon_positions

def find_splice_sites(seq, start, end, orientation):
    """Adjust the intron boundaries based on the closest splice site (CT-AC for antisense)."""
    SEARCH_RANGE = 20  # Define the search range for splice sites

    if orientation == '+':
        # For sense strands (GT-AG logic)
        new_start = seq.find('GT', max(0, start - SEARCH_RANGE), start)  # Search for GT donor site
        new_end = seq.rfind('AG', end, min(len(seq), end + SEARCH_RANGE))  # Search for AG acceptor site
    else:
        # For antisense strands (CT-AC logic)
        new_start = seq.find('CT', max(0, start - SEARCH_RANGE), start)  # Search for CT donor site
        new_end = seq.rfind('AC', end, min(len(seq), end + SEARCH_RANGE))  # Search for AC acceptor site

    # Adjust start and end if valid splice sites are found
    if new_start != -1:
        start = new_start + 1  # Adjust for 1-based indexing
    if new_end != -1:
        end = new_end + 2  # Adjust for length of AG or AC

    return start, end

def check_start_stop_codons(exon_positions, seq, orientation):
    """Check for start and stop codons in the exons."""
    # Debugging: print exon_positions and sequence
    print(f"Exon positions: {exon_positions}")
    print(f"Sequence length: {len(seq)}")

    if orientation == '+':
        start_codon = seq[exon_positions[0][1] - 1:exon_positions[0][1] + 2]
        
        if start_codon == 'ATG':
            print("Start codon (ATG) found at the start of exon 1.")
        else:
            print(f"No start codon found. Sequence start: {start_codon}")
        
        stop_codon = seq[exon_positions[-1][2]:exon_positions[-1][2] + 3]

        if stop_codon in ['TAA', 'TAG', 'TGA']:
            print(f"Stop codon ({stop_codon}) found at the end of exon 3.")
        else:
            print(f"No stop codon found.")

    else:
        # Antisense strand
        if exon_positions[-1][2] is not None:
            start_codon = seq[exon_positions[-1][2] - 3:exon_positions[-1][2]]
            
            if start_codon == 'TAC':
                print("Start codon (TAC) found at the end of exon 3.")
            else:
                print(f"No start codon found. Sequence start: {start_codon}")
        else:
            print("Error: exon_positions[-1][2] is None.")

        if exon_positions[0][1] is not None:
            stop_codon = seq[exon_positions[0][1] - 4:exon_positions[0][1] - 1]
            print(f"Stop codon found: {stop_codon}")
            if stop_codon in ['TTA', 'CTA', 'TCA']:
                print(f"Stop codon ({stop_codon}) found before exon 1.")
            else:
                print(f"No stop codon found.")
        else:
            print("Error: exon_positions[0][1] is None.")


def write_gff_file(gff_file, exon_positions, intron_positions, sequence_id, protein_id):
    """Write the exon and intron positions to a new GFF file."""
    with open(gff_file, "w") as gff:
        for exon_name, start, end, orientation in exon_positions:
            gff.write(f"{sequence_id}\tgenbank\t{exon_name}\t{start}\t{end}\t.\t{orientation}\t.\tProtein ID: {protein_id}; exon_number={exon_name[-1]}\n")
        
        for intron_name, start, end, orientation in intron_positions:
            gff.write(f"{sequence_id}\tgenbank\t{intron_name}\t{start}\t{end}\t.\t.\t.\tProtein ID: {protein_id}; intron_number={intron_name[-1]}\n")


def write_introns_to_fasta(intron_positions, exon_positions, seq, fasta_prefix="intron"):
    """Write the intron sequences to separate FASTA files based on their positions, without reverse complementing."""
    
    # Check if orientation is antisense
    if exon_positions[0][3] == "-":
        # Reverse the order of the exons and introns
        exon_positions = exon_positions[::-1]
        intron_positions = intron_positions[::-1]
    
    for idx, (intron_name, intron_start, intron_end, orientation) in enumerate(intron_positions, start=1):
        # Extract the intron sequence from the full sequence
        intron_seq = seq[intron_start-1:intron_end]  # 1-based to 0-based index adjustment

        # Write the sequence without reverse complementing for antisense
        fasta_filename = f"{fasta_prefix}{idx}.fasta"
        with open(fasta_filename, "w") as fasta_file:
            fasta_file.write(f">{intron_name} {intron_start}-{intron_end} {orientation}\n")
            fasta_file.write(format_sequence(intron_seq))

        print(f"Intron {idx} written to {fasta_filename}")

def write_exons_to_fasta(exon_positions, seq, fasta_prefix="exon"):
    """Write the exon sequences to separate FASTA files."""
    
    for idx, (exon_name, exon_start, exon_end, orientation) in enumerate(exon_positions, start=1):
        # Extract the exon sequence from the full sequence
        exon_seq = seq[exon_start-1:exon_end]  # 1-based to 0-based index adjustment

        # Write the sequence without reverse complementing for antisense
        fasta_filename = f"{fasta_prefix}{idx}.fasta"
        with open(fasta_filename, "w") as fasta_file:
            fasta_file.write(f">{exon_name} {exon_start}-{exon_end} {orientation}\n")
            fasta_file.write(format_sequence(exon_seq))

        print(f"Exon {idx} written to {fasta_filename}")

def format_sequence(sequence, line_length=60):
    """Format the sequence to have a specific number of characters per line (default 60)."""
    return '\n'.join(sequence[i:i + line_length] for i in range(0, len(sequence), line_length))

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
    intron_positions = adjust_intron_boundaries(exon_positions, full_seq)
    exon_positions = adjust_exon_boundaries(exon_positions, intron_positions, exon_positions[0][3])
    print(f"Exon positions after adjustment: {exon_positions}")  # Debugging output

    # Step 5: Check for start and stop codons
    if exon_positions is not None:
        check_start_stop_codons(exon_positions, full_seq, orientation)
    else:
        print("Error: exon_positions is None. Cannot check start/stop codons.")

    # Step 6: Write new GFF file
    output_gff_file = f"{protein_id}_chk.gff"
    write_gff_file(output_gff_file, exon_positions, intron_positions, sequence_id, protein_id)

    # Step 7: Write exon and intron FASTA files
    write_introns_to_fasta(intron_positions, exon_positions, full_seq)
    write_exons_to_fasta(exon_positions, full_seq)
    #excise_exons_introns_to_fasta(exon_positions, intron_positions, full_seq)

    print(f"Adjusted GFF file '{output_gff_file}' and exon/intron FASTA files have been created successfully.")

if __name__ == "__main__":
    main()
