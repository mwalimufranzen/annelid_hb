def find_start_codon(nucleotide_seq, orientation, exon_positions):
    """Find the correct start codon (ATG for sense or CAT for antisense).
       If not found, search for the start codon and adjust exon positions."""
    
    if orientation == "+":
        # Check the first exon for "ATG"
        first_exon_start = exon_positions[0][1] - 1  # 0-based index
        if nucleotide_seq[first_exon_start:first_exon_start + 3] != "ATG":  # If the first exon does not start with ATG
            print("Start codon (ATG) not found at the beginning of the first exon. Searching upstream...")
            # Search upstream from the first exon
            upstream_region = nucleotide_seq[max(0, first_exon_start - 100):first_exon_start]  # Search up to 100 bases upstream
            start_index = upstream_region.rfind("ATG")  # Look for the last ATG before exon 1

            if start_index != -1:
                added_bases = first_exon_start - (first_exon_start - 100 + start_index)  # Calculate the number of added bases
                print(f"Start codon (ATG) found {added_bases} bases upstream. Adjusting exon 1.")
                # Adjust exon 1 position
                exon_positions[0] = (exon_positions[0][0], exon_positions[0][1] - added_bases, exon_positions[0][2], exon_positions[0][3])

                # Update the exon1.fasta file
                exon_seq = nucleotide_seq[first_exon_start - added_bases:exon_positions[0][2]]
                with open(f"exon1_sense.fasta", "w") as exon_file:
                    exon_file.write(f">exon1 {exon_positions[0][1]}-{exon_positions[0][2]} +\n")
                    exon_file.write(format_sequence(exon_seq))
                print("exon1.fasta file has been updated.")
                
                return nucleotide_seq[added_bases:], exon_positions  # Return adjusted sequence
            else:
                print("No start codon (ATG) found in the upstream region. Exiting...")
                sys.exit(1)
        else:
            print("Start codon (ATG) found at the beginning of the first exon.")
    
    elif orientation == "-":
        # For antisense, check exon 1 for "CAT" at the end
        first_exon_end = exon_positions[0][2]
        if nucleotide_seq[first_exon_end - 3:first_exon_end] != "CAT":  # If the first exon does not end with CAT
            print("Stop codon (CAT) not found at the end of the first exon. Searching downstream...")
            # Search downstream after the first exon
            downstream_region = nucleotide_seq[first_exon_end:min(len(nucleotide_seq), first_exon_end + 100)]  # Search up to 100 bases downstream
            end_index = downstream_region.find("CAT")

            if end_index != -1:
                print(f"Stop codon (CAT) found {end_index + 1} bases downstream. Adjusting exon 1.")
                added_bases = end_index + 1
                # Adjust the first exon position
                exon_positions[0] = (exon_positions[0][0], exon_positions[0][1], exon_positions[0][2] + added_bases, exon_positions[0][3])
                
                # Update the first exon fasta file for antisense and sense
                exon_seq = nucleotide_seq[exon_positions[0][1] - 1:exon_positions[0][2]]
                
                # Output the original sequence as exon#_sense
                with open(f"exon1_sense.fasta", "w") as exon_file:
                    exon_file.write(f">exon1 {exon_positions[0][1]}-{exon_positions[0][2]} +\n")
                    exon_file.write(format_sequence(exon_seq))
                print(f"exon1_sense.fasta has been updated.")
                
                # Output the reverse complement as exon#_antisense
                reverse_complement_seq = Seq(exon_seq).reverse_complement()
                with open(f"exon1_antisense.fasta", "w") as exon_file:
                    exon_file.write(f">exon1 {exon_positions[0][1]}-{exon_positions[0][2]} -\n")
                    exon_file.write(format_sequence(str(reverse_complement_seq)))
                print(f"exon1_antisense.fasta has been updated.")
                
                return nucleotide_seq[:first_exon_end + end_index], exon_positions  # Return adjusted sequence
            else:
                print("No stop codon (CAT) found in the downstream region. Exiting...")
                sys.exit(1)
        else:
            print("Stop codon (CAT) found at the end of the first exon.")

    return nucleotide_seq, exon_positions  # Return original if no changes