from Bio.Seq import Seq
from Bio import SeqIO

def translate_fasta(fasta_file):
    print("The script accepts input in a file named input_cDNA.fasta.")
    
    # Open the input FASTA file
    with open(fasta_file, "r") as input_handle:
        for record in SeqIO.parse(input_handle, "fasta"):
            nucleotide_seq = record.seq
            seq_length = len(nucleotide_seq)
            
            # Count the number of nucleotides and check if it's a multiple of 3
            is_multiple_of_three = (seq_length % 3 == 0)
            
            print(f"Processing sequence {record.id}:")
            print(f"Nucleotide count: {seq_length}")
            print(f"Is multiple of 3: {'Yes' if is_multiple_of_three else 'No'}")
            
            # Translate three possible open reading frames (ORFs) from the forward strand
            for frame in range(3):
                translated_seq = nucleotide_seq[frame:].translate()

                # Create an ORF-specific FASTA file for forward strand
                orf_filename = f"{record.id}_ORF{frame + 1}_forward.fasta"
                with open(orf_filename, "w") as orf_handle:
                    orf_handle.write(f">{record.id}_ORF{frame + 1}_forward\n")
                    orf_handle.write(f"{translated_seq}\n")
                
                print(f"Open Reading Frame {frame + 1} (Forward):")
                print(translated_seq)

            # Get the reverse complement of the sequence
            reverse_complement_seq = nucleotide_seq.reverse_complement()

            # Translate three possible open reading frames (ORFs) from the reverse complement strand
            for frame in range(3):
                translated_seq = reverse_complement_seq[frame:].translate()

                # Create an ORF-specific FASTA file for reverse strand
                orf_filename = f"{record.id}_ORF{frame + 1}_reverse.fasta"
                with open(orf_filename, "w") as orf_handle:
                    orf_handle.write(f">{record.id}_ORF{frame + 1}_reverse\n")
                    orf_handle.write(f"{translated_seq}\n")
                
                print(f"Open Reading Frame {frame + 1} (Reverse):")
                print(translated_seq)

# Example usage
fasta_file = "input_cDNA.fasta"  # Replace with your input FASTA file
translate_fasta(fasta_file)
