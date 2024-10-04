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
    match = re.search(r'Query #1: tr\|(\w+)\|', content)
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
        # Extract only the first Sbjct line after each range using a more flexible regex
        sbjct_match = re.search(r'Sbjct\s+(\d+)\s+[^\n]+\s+(\d+)', content[match.end():])
        if sbjct_match:
            sbjct_start, sbjct_end = map(int, sbjct_match.groups())
            # Determine strand based on the first Sbjct number
            if sbjct_start == start:
                orientation = '+'
            elif sbjct_start == end:
                orientation = '-'
            else:
                sys.stderr.write(f"Unable to determine strand orientation for range {start}-{end}.\n")
                sys.stderr.write(f"Sbjct start: {sbjct_start}, Sbjct end: {sbjct_end}, expected start: {start}, expected end: {end}\n")
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


def calculate_intron_positions(exon_positions):
    """Calculate intron positions from exon positions."""
    intron_positions = []
    exon_list = [(start, end) for _, start, end, _ in exon_positions]

    if len(exon_list) >= 2:
        intron_positions.append(('intron1', exon_list[0][1] + 1, exon_list[1][0] - 1))
    if len(exon_list) == 3:
        intron_positions.append(('intron2', exon_list[1][1] + 1, exon_list[2][0] - 1))

    return intron_positions


def write_gff_file(gff_file, exon_positions, intron_positions, sequence_id, protein_id):
    """Write the exon and intron positions to a GFF file."""
    with open(gff_file, "w") as gff:
        # Write exons
        for exon_name, start, end, orientation in exon_positions:
            gff.write(f"{sequence_id}\tgenbank\t{exon_name}\t{start}\t{end}\t.\t{orientation}\t0\tProtein ID: {protein_id}; exon_number={exon_name[-1]}\n")
        
        # Write introns
        for intron_name, start, end in intron_positions:
            gff.write(f"{sequence_id}\tgenbank\t{intron_name}\t{start}\t{end}\t.\t.\t0\tProtein ID: {protein_id}; intron_number={intron_name[-1]}\n")


def excise_exons_to_fasta(exon_positions, target_seq, protein_id):
    """Excise the exons from the genomic data and write them to exon#_sense.fasta files."""
    for exon_name, start, end, orientation in exon_positions:
        exon_seq = target_seq[start - 1:end]  # Adjust for 0-based indexing
        # Write exon to file based on orientation
        with open(f"{exon_name}_sense.fasta", "w") as exon_file:
            exon_file.write(f">{exon_name}_sense {start}-{end}\n{exon_seq}\n")


def write_intron_files(exon_positions, seq):
    """Extract and write introns based on exon positions, with sense/antisense designations."""
    if len(exon_positions) > 1:
        for j in range(1, len(exon_positions)):
            intron_start = exon_positions[j-1][2] + 1  # exon end of previous exon
            intron_end = exon_positions[j][1] - 1  # exon start of current exon
            if intron_start < intron_end:
                intron_seq = seq[intron_start-1:intron_end]
                orientation = exon_positions[j-1][3]  # Use the orientation of the previous exon
                intron_filename = f"intron{j}_{'sense' if orientation == '+' else 'antisense'}.fasta"
                with open(intron_filename, 'w') as intron_file:
                    SeqIO.write(SeqRecord(Seq(intron_seq), id=f"intron{j}", description=f"intron{j} range {intron_start}-{intron_end}"), intron_file, 'fasta')
                sys.stdout.write(f"Intron {j} written to {intron_filename}.\n")


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

    # Calculate intron positions based on exon positions
    intron_positions = calculate_intron_positions(exon_positions)

    # Generate the GFF file using the Protein ID for naming
    gff_file = f"{protein_id}.gff"
    write_gff_file(gff_file, exon_positions, intron_positions, sequence_id, protein_id)

    # Excise the exons from the target sequence and write to files
    excise_exons_to_fasta(exon_positions, target_seq, protein_id)

    # Write intron files
    write_intron_files(exon_positions, target_seq)

    print(f"GFF file '{gff_file}', target sequence, and intron/exon FASTA files have been created successfully.")


if __name__ == "__main__":
    main()
