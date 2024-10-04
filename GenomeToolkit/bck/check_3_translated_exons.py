from Bio import SeqIO
from Bio.Seq import Seq

# File names for exons
exon_files = ['exon1.fasta', 'exon2.fasta', 'exon3.fasta']

# Output file names
combined_output = 'open_reading_frame.fasta'
translated_output = 'translated_reading_frame.fasta'

# Combine exon sequences into a single contiguous sequence
combined_sequence = ""

for exon_file in exon_files:
    # Reading the exon sequences
    with open(exon_file, 'r') as handle:
        for record in SeqIO.parse(handle, "fasta"):
            combined_sequence += str(record.seq)

# Write combined sequence to a new FASTA file
with open(combined_output, 'w') as combined_handle:
    combined_handle.write(">open_reading_frame\n")
    combined_handle.write(combined_sequence + "\n")

# Translate the combined sequence (first open reading frame)
combined_seq_obj = Seq(combined_sequence)
translated_sequence = combined_seq_obj.translate(to_stop=True)  # Stop at first stop codon

# Write translated sequence to a new FASTA file
with open(translated_output, 'w') as translated_handle:
    translated_handle.write(">translated_reading_frame\n")
    translated_handle.write(str(translated_sequence) + "\n")

print(f"Combined sequence saved to {combined_output}")
print(f"Translated sequence saved to {translated_output}")
