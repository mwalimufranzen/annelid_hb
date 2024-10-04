from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def reverse_transcribe_fasta(input_fasta, reverse_transcribed_fasta, translated_fasta):
    # Read the nucleotide sequences from the input FASTA file
    records = list(SeqIO.parse(input_fasta, "fasta"))

    # Prepare lists to store reverse transcribed and translated sequences
    reverse_transcribed_records = []
    translated_records = []

    for record in records:
        # Reverse transcribe the sequence (RNA -> cDNA)
        reverse_transcribed_seq = record.seq.reverse_complement()
        reverse_transcribed_record = SeqRecord(
            reverse_transcribed_seq,
            id=record.id,
            description="Reverse transcribed cDNA sequence"
        )
        reverse_transcribed_records.append(reverse_transcribed_record)

        # Translate the cDNA sequence to amino acids
        translated_seq = reverse_transcribed_seq.translate(to_stop=True)
        translated_record = SeqRecord(
            translated_seq,
            id=record.id,
            description="Translated amino acid sequence"
        )
        translated_records.append(translated_record)

    # Write the reverse transcribed cDNA sequences to the output FASTA file
    SeqIO.write(reverse_transcribed_records, reverse_transcribed_fasta, "fasta")

    # Write the translated amino acid sequences to the output FASTA file
    SeqIO.write(translated_records, translated_fasta, "fasta")

if __name__ == "__main__":
    input_fasta = "input_nucleotide.fasta"
    reverse_transcribed_fasta = "reverse_transcribed.fasta"
    translated_fasta = "translated.fasta"
    
    reverse_transcribe_fasta(input_fasta, reverse_transcribed_fasta, translated_fasta)
