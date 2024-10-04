def reverse_complement(seq):
    # Create a dictionary to map each nucleotide to its complement
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 
                  'a': 't', 't': 'a', 'c': 'g', 'g': 'c'}  # Lowercase support

    # Replace each nucleotide with its complement
    comp_seq = ''.join([complement[base] for base in seq])

    # Return the reverse of the complemented sequence
    return comp_seq[::-1]

# Example usage
sequence = "CCAGACCGTTACTTCGA"
rev_comp = reverse_complement(sequence)
print(f"Original sequence: {sequence}")
print(f"Reverse complement: {rev_comp}")
