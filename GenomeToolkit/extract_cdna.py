import os
import requests
from Bio import SeqIO

def download_fasta_from_ncbi(sequence_id, local_file_path):
    """
    Downloads a FASTA file from NCBI and saves it locally.
    """
    # Define the URL for downloading the FASTA file
    url = f"https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id={sequence_id}&db=nucleotide&report=fasta&retmode=text"

    # Send a request to fetch the FASTA data
    response = requests.get(url)
    if response.status_code == 200:
        with open(local_file_path, 'w') as out_file:
            out_file.write(response.text)
    else:
        raise RuntimeError(f"Failed to download FASTA file. HTTP status code: {response.status_code}")

def extract_cdna_from_alignment(alignment_file, output_fasta_file):
    """
    Extracts a cDNA sequence from a transcriptome fragment based on the alignment file
    and writes the sequence to a FASTA file.
    """
    # Read alignment file to get the sequence ID and range limits
    with open(alignment_file, 'r') as f:
        alignment_lines = f.readlines()

    # Parse the alignment file
    sequence_id = None
    start = None
    end = None
    
    for line in alignment_lines:
        if line.startswith('Sequence ID:'):
            sequence_id = line.split(':')[1].strip()
        elif line.startswith('Range 1:'):
            range_part = line.split(':')[1].strip()
            start, end = map(int, range_part.split(' to '))
    
    if not sequence_id or start is None or end is None:
        raise ValueError("Unable to parse alignment file correctly.")

    # Define the local file path to save the FASTA file
    fasta_file_path = f"{sequence_id}.fasta"
    
    # Download the FASTA file from NCBI
    download_fasta_from_ncbi(sequence_id, fasta_file_path)
    
    # Read the downloaded FASTA file
    with open(fasta_file_path, 'r') as f:
        record = SeqIO.read(f, "fasta")
    
    # Extract the cDNA sequence using the specified range
    cDNA_sequence = record.seq[start-1:end]  # start-1 for zero-based index

    # Write the cDNA sequence to a FASTA file
    with open(output_fasta_file, 'w') as f:
        f.write(f">{sequence_id}\n")
        for i in range(0, len(cDNA_sequence), 80):  # Optional: wrap lines at 80 characters
            f.write(str(cDNA_sequence[i:i+80]) + '\n')  # Convert Seq to str before writing
    
    # Optionally remove the downloaded FASTA file
    os.remove(fasta_file_path)

if __name__ == "__main__":
    # Specify the filenames
    alignment_file = 'alignment.txt'   # File containing alignment details
    output_fasta_file = 'extracted_cDNA.fasta'  # Output FASTA file
    
    # Extract the cDNA sequence and write to FASTA
    extract_cdna_from_alignment(alignment_file, output_fasta_file)
