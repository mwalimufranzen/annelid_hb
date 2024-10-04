from Bio.Seq import Seq

def translate_sequence(nucleotide_seq, orientation):
    """Translate the sequence based on the orientation (sense or antisense)"""
    nucleotide_seq = Seq(nucleotide_seq)  # Ensure nucleotide_seq is a Seq object
    
    # Check orientation and handle translation
    if orientation == "+":
        orf = nucleotide_seq  # Use the sequence as is for sense strand
    elif orientation == "-":
        orf = nucleotide_seq.reverse_complement()  # Use the reverse complement for antisense strand
    else:
        raise ValueError("Invalid orientation: expected '+' or '-'")
    
    # Translate to amino acid sequence, stopping at the first stop codon
    amino_acid_seq = orf.translate(to_stop=True)
    
    return str(amino_acid_seq)

# Example usage with the correct nucleotide sequence and orientation
nucleotide_seq = """
GCCTTCCAGAGGACGTGTCCGAACTCGTCCCTGACCCTGCTGGTGCCGTACACCTTGTTC
CACTGGGTCTTGACCAGGAGGCGCTGCATGGGGCCGCACTCTGCAGCCACCATGGCCACC
AGGGCGAACAGAACCACCAAGAACTTCATATCGAAGTAGCGGGTTGGGATGTGCCTTTCT
ACGTGCTGCACGTGAAGGTGAGCCAGGGCGACGTCAAGGTCAGCCTGGTTGTCAAGGAGG
GAGATGGCAATGTCAAGGCCACCCAAGACACGGGCGCTGTGAGCCATGAACTCGGGAGAG
TAGATGTCGTCACCGTTGACTCTCTTGAAGAGAGCCCGGGTCTCGGGATCCTGGGCGAAG
ACACTTGATGCCGTTGGCGATGATGTCAAAGCACGAGCTCCAGGCGGTCTTATCGAAGCA
GCGTCCCAGGGCGCTGGGGGCATACTCCATCAGGGCGTTCTTGAACAGCTATAGAGATAA
AATATACACGGAATGAATGAATCGGAGGTTTAGGCAGTTTTGTGGAAACAAGGGCTAAGT
CGCGCCAAGAAACAT""".replace('\n', '')

# Test for antisense translation
print(translate_sequence(nucleotide_seq, "-"))
