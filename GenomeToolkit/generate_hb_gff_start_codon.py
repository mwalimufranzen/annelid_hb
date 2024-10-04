import re
import sys
import requests
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def download_sequence(sequence_id, output_file):
    """Download the sequence from NCBI using the sequence ID and save it to an output file."""
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

def parse_alignment_file(alignment_file):
    """Parse the alignment.txt file to extract relevant information."""
    exon_positions = []
    sequence_id = None
    protein_id = None
    content = None

    # Read the alignment file content
    with open(alignment_file, 'r') as f:
        content = f.read()

    # Extract Protein ID
    match = re.search(r'Protein ID: (\w+)', content)
    if match:
        protein_id = match.group(1)
        print(f"Protein ID found: {protein_id}")

    # Extract Sequence ID (without "Length")
    match = re.search(r'Sequence ID: (\w+\.\d+)', content)
    if match:
        sequence_id = match.group(1)
        print(f"Sequence ID found: {sequence_id}")

    # Extract ranges and determine strand orientation based on the first occurrence of Sbjct line
    ranges = []
    for match in re.finditer(r'Range \d+: (\d+) to (\d+)', content):
        start, end = map(int, match.groups())
        sbjct_match = re.search(r'Sbjct\s+(\d+)\s+[^\n]+\s+(\d+)', content[match.end():])
        if sbjct_match:
            sbjct_start, sbjct_end = map(int, sbjct_match.groups())
            if sbjct_start == start:
                orientation = '+'
            elif sbjct_start == end:
                orientation = '-'
            else:
                sys.stderr.write(f"Unable to determine strand orientation for range {start}-{end}.\n")
                sys.exit(1)
            ranges.append((start, end, orientation))
            sys.stdout.write(f"Range {len(ranges)}: start {start}, end {end}, {orientation}\n")
        else:
            sys.stderr.write(f"No Sbjct line found for range {start}-{end}. Trying to proceed.\n")
            ranges.append((start, end, 'unknown'))  # Allow it to proceed even without orientation

    # Store exon positions
    for idx, (start, end, orientation) in enumerate(ranges, 1):
        exon_name = f"exon{idx}"
        exon_positions.append((exon_name, start, end, orientation))
    return protein_id, sequence_id, exon_positions

def find_start_codon(nucleotide_seq, orientation, exon_positions):
    """Find the correct start codon (ATG for sense or CAT for antisense).
       If not found, search for the start codon and adjust exon positions."""
    
    if orientation == "+":
        # Check the first exon for "ATG"
        first_exon_start = exon_positions[0][1] - 1  # 0-based index
        if nucleotide_seq[first_exon_start:first_exon_start + 3] != "ATG":  # If the first exon does not start with ATG
            print("Start codon (ATG) not found at the beginning of the first exon. Searching upstream...")
            # Search upstream from the first exon (up to 60 bases upstream)
            upstream_region = nucleotide_seq[max(0, first_exon_start - 60):first_exon_start]
            start_index_upstream = upstream_region.rfind("ATG")
            
            if start_index_upstream != -1:
                # Start codon found upstream, adjust exon
                added_bases = first_exon_start - (first_exon_start - 60 + start_index_upstream)
                print(f"Start codon (ATG) found {added_bases} bases upstream. Adjusting exon 1.")
                exon_positions[0] = (exon_positions[0][0], exon_positions[0][1] - added_bases, exon_positions[0][2], exon_positions[0][3])

            else:
                # Search downstream (up to 40 bases downstream)
                print("Searching downstream for a closer start codon...")
                downstream_region = nucleotide_seq[first_exon_start:first_exon_start + 40]
                start_index_downstream = downstream_region.find("ATG")
                
                if start_index_downstream != -1:
                    # Start codon found downstream, adjust exon
                    exon_positions[0] = (exon_positions[0][0], exon_positions[0][1] + start_index_downstream, exon_positions[0][2], exon_positions[0][3])
                    print(f"Start codon (ATG) found {start_index_downstream + 1} bases downstream. Shortening exon 1.")
                else:
                    print("Error: No valid start codon (ATG) found within the allowed range. Exiting...")
                    sys.exit(1)

    elif orientation == "-":
        # For antisense, check exon 1 for "CAT" at the end
        first_exon_end = exon_positions[0][2]
        if nucleotide_seq[first_exon_end - 3:first_exon_end] != "CAT":  # If the first exon does not end with CAT
            print("Start codon (CAT) not found at the end of the first exon. Searching downstream...")
            # Search downstream (up to 60 bases downstream from the highest position)
            downstream_region = nucleotide_seq[first_exon_end:first_exon_end + 60]
            stop_index_downstream = downstream_region.find("CAT")
            
            if stop_index_downstream != -1:
                # Stop codon found downstream, adjust exon
                print(f"Start codon (CAT) found at position {first_exon_end + stop_index_downstream}.")
                exon_positions[0] = (exon_positions[0][0], exon_positions[0][1], first_exon_end + stop_index_downstream + 3, exon_positions[0][3])
            else:
                # Search upstream (up to 40 bases upstream)
                print("Searching upstream for a closer stop codon...")
                upstream_region = nucleotide_seq[max(0, first_exon_end - 40):first_exon_end]
                stop_index_upstream = upstream_region.rfind("CAT")
                
                if stop_index_upstream != -1:
                    # Stop codon found upstream, adjust exon
                    print(f"Start codon (CAT) found upstream at position {first_exon_end - (40 - stop_index_upstream)}.")
                    exon_positions[0] = (exon_positions[0][0], exon_positions[0][1], first_exon_end - (40 - stop_index_upstream), exon_positions[0][3])
                else:
                    print("Error: No valid start codon (CAT) found within the allowed range. Exiting...")
                    sys.exit(1)

    # Write the corrected exon sequences
    write_exon_fasta(nucleotide_seq, exon_positions[0][1], exon_positions[0][2], orientation, 1)
    
    return nucleotide_seq, exon_positions

def write_exon_fasta(nucleotide_seq, start, end, orientation, exon_number):
    """Write the exon sequence to a FASTA file with correct orientation."""
    exon_seq = nucleotide_seq[start - 1:end]  # Adjust for 0-based indexing
    
    # Write both sense and antisense sequences to separate files
    with open(f"exon{exon_number}_sense.fasta", "w") as exon_file:
        exon_file.write(f">exon{exon_number} {start}-{end} +\n")
        exon_file.write(format_sequence(exon_seq))
    
    if orientation == "-":
        # Write the reverse complement for the antisense case
        reverse_complement_seq = Seq(exon_seq).reverse_complement()
        with open(f"exon{exon_number}_antisense.fasta", "w") as exon_file:
            exon_file.write(f">exon{exon_number} {start}-{end} -\n")
            exon_file.write(format_sequence(str(reverse_complement_seq)))

def format_sequence(sequence):
    """Helper function to format sequence into 60 characters per line for FASTA format."""
    return '\n'.join(sequence[i:i + 60] for i in range(0, len(sequence), 60))

def excise_exons_introns_to_fasta(exon_positions, target_seq):
    """Excise exons and introns from the genomic data and write them to exon#.fasta and intron#.fasta files.
       Outputs in standard FASTA format (60 nucleotides per line) with sequence range in the header."""
    
    def format_sequence(sequence):
        """Format the sequence to 60 characters per line."""
        return '\n'.join(sequence[i:i+60] for i in range(0, len(sequence), 60))

    # Write exons
    for idx, (exon_name, start, end, orientation) in enumerate(exon_positions):
        exon_seq = target_seq[start - 1:end]  # Adjust for 0-based indexing
        
        #if orientation == '+':
        exon_header = f">exon{idx + 1} {start}-{end} {orientation}"  # Include range in the header
        formatted_exon_seq = format_sequence(exon_seq)

        with open(f"exon{idx + 1}_sense.fasta", "w") as exon_file:
            exon_file.write(f"{exon_header}\n{formatted_exon_seq}\n")
        print(f"Exon {idx + 1} written: start={start}, end={end}")
        
        #elif orientation == '-':
        if orientation == '-':
            # Reverse complement the exon sequence for antisense strand
            reverse_complement_seq = Seq(exon_seq).reverse_complement()
            antisense_header = f">exon{idx + 1}_antisense {start}-{end} {orientation}"
            formatted_antisense_seq = format_sequence(str(reverse_complement_seq))

            with open(f"exon{idx + 1}_antisense.fasta", "w") as exon_file:
                exon_file.write(f"{antisense_header}\n{formatted_antisense_seq}\n")
            #print(f"Reverse complement exon {idx + 1} written.")

    # Write introns
    for idx in range(1, len(exon_positions)):
        exon_list = [(start, end, orientation) for _, start, end, orientation in exon_positions]

        if exon_list[idx][2] == '-':  # Handle antisense strand (reverse intron order)
            intron_start = exon_list[idx][1] + 1  # End of next exon + 1
            intron_end = exon_list[idx-1][0] - 1  # Start of current exon - 1
            intron_seq = target_seq[intron_start - 1:intron_end]  # Adjust for 0-based indexing
            reverse_complement_intron = Seq(intron_seq).reverse_complement()  # Reverse complement the intron

            intron_header = f">intron{idx}_antisense {intron_start}-{intron_end}"
            formatted_intron_seq = format_sequence(str(reverse_complement_intron))

            with open(f"intron{idx}_antisense.fasta", "w") as intron_file:
                intron_file.write(f"{intron_header}\n{formatted_intron_seq}\n")
            #print(f"Reverse complement intron {idx} written: start={intron_start}, end={intron_end}")

        #else:  # Handle sense strand
        if exon_list[idx][2] == '+':
            intron_start = exon_list[idx-1][1] + 1  # End of previous exon + 1
            intron_end = exon_list[idx][0] - 1  # Start of next exon - 1
        else:
            intron_start = exon_list[idx][1] + 1  # End of next exon + 1
            intron_end = exon_list[idx-1][0] - 1  # Start of current exon - 1

        if intron_start < intron_end:
            intron_seq = target_seq[intron_start - 1:intron_end]  # Adjust for 0-based indexing
            intron_header = f">intron{idx}_sense {intron_start}-{intron_end}"
            formatted_intron_seq = format_sequence(intron_seq)

            with open(f"intron{idx}_sense.fasta", "w") as intron_file:
                intron_file.write(f"{intron_header}\n{formatted_intron_seq}\n")
            print(f"Intron {idx} written: start={intron_start}, end={intron_end}")
        else:
            print(f"Warning: intron_start ({intron_start}) > intron_end ({intron_end}) for intron {idx}")

def calculate_intron_positions(exon_positions):
    """Calculate intron positions from exon positions based on the specified scheme."""
    intron_positions = []
    exon_list = [(start, end, orientation) for _, start, end, orientation in exon_positions]

    # Iterate through exons to calculate introns between them
    for i in range(1, len(exon_list)):
        intron_start = exon_list[i-1][1] + 1  # End of previous exon + 1
        intron_end = exon_list[i][0] - 1  # Start of current exon - 1
        intron_name = f"intron{i}"

        # Ensure valid intron range
        if intron_start <= intron_end:
            intron_positions.append((intron_name, intron_start, intron_end, exon_list[i-1][2]))

    return intron_positions


def write_gff_file(gff_file, exon_positions, intron_positions, sequence_id, protein_id):
    """Write the exon and intron positions to a GFF file with correct start and end."""
    with open(gff_file, "w") as gff:
        # Write exons
        for exon_name, start, end, orientation in exon_positions:
            gff.write(f"{sequence_id}\tgenbank\t{exon_name}\t{min(start, end)}\t{max(start, end)}\t.\t{orientation}\t0\tProtein ID: {protein_id}; exon_number={exon_name[-1]}\n")
        
        # Write introns between exons
        for intron_name, start, end, orientation in intron_positions:  # Added 'orientation' to match the tuple
            gff.write(f"{sequence_id}\tgenbank\t{intron_name}\t{min(start, end)}\t{max(start, end)}\t.\t{orientation}\t0\tProtein ID: {protein_id}; intron_number={intron_name[-1]}\n")

def create_target_with_flanking(seq, exon_positions):
    """Create a target sequence with 350 bp flanking regions."""
    min_pos = min([start for _, start, _, _ in exon_positions])
    max_pos = max([end for _, _, end, _ in exon_positions])

    # Add 350 bp flanking sequences
    flanking_start = max(min_pos - 350, 0)
    flanking_end = min(max_pos + 350, len(seq))
    target_seq = seq[flanking_start:flanking_end]

    target_filename = "target.fasta"
    with open(target_filename, 'w') as file:
        target_record = SeqRecord(Seq(target_seq),
                                  id="target",
                                  description="Target sequence with 350bp flanking regions")
        SeqIO.write(target_record, file, 'fasta')

    sys.stdout.write(f"Target sequence written to {target_filename}.\n")
    
    # Adjust exon positions relative to the new target sequence
    for i in range(len(exon_positions)):
        exon_positions[i] = (exon_positions[i][0], exon_positions[i][1] - flanking_start, exon_positions[i][2] - flanking_start, exon_positions[i][3])
    
    return target_seq, exon_positions


def main():
    alignment_file = "alignment.txt"
    
    # Parse the alignment file to get Protein ID, Sequence ID, and exon positions
    protein_id, sequence_id, exon_positions = parse_alignment_file(alignment_file)
    
    # Check if any of the necessary information is missing and provide detailed error messages
    if not protein_id:
        print("Error: Protein ID not found in the alignment file.")
    if not sequence_id:
        print("Error: Sequence ID not found in the alignment file.")
    if not exon_positions:
        print("Error: Exon positions not found in the alignment file.")
        return

    # Download the sequence from NCBI using the Sequence ID
    fasta_file = f"{sequence_id}.fasta"
    download_sequence(sequence_id, fasta_file)

    # Read the downloaded sequence
    target_record = SeqIO.read(fasta_file, "fasta")
    full_seq = str(target_record.seq)

    # Create the target sequence with 350 bp flanking regions and adjust exon positions
    target_seq, exon_positions = create_target_with_flanking(full_seq, exon_positions)

    # ** Find start codon (ATG for sense, CAT for antisense) and adjust exon positions if needed **
    orientation = exon_positions[0][3]  # Get orientation from the first exon
    # Call find_start_codon after creating the target sequence and reading exon positions
    target_seq, exon_positions = find_start_codon(target_seq, orientation, exon_positions)

    # Calculate intron positions based on exon positions
    intron_positions = calculate_intron_positions(exon_positions)

    # Generate the GFF file using the Protein ID for naming
    gff_file = f"{protein_id}.gff"
    write_gff_file(gff_file, exon_positions, intron_positions, sequence_id, protein_id)

    # Excise the exons from the target sequence and write to files
    excise_exons_introns_to_fasta(exon_positions, target_seq)

    print(f"GFF file '{gff_file}', target sequence, and intron/exon FASTA files have been created successfully.")

if __name__ == "__main__":
    main()
