from re import findall
from time import perf_counter

# Generic, for loop solution
def count_kmer_loop(sequence, kmer):

    kmer_count = 0
    for position in range(len(sequence) - 1):
        if sequence[position:position + len(kmer)] == kmer:
            kmer_count += 1
    return kmer_count
"""
if sequence[position:position + len(kmer)] == kmer

if ATT[TAA]AAC == AAA | do nothing
if ATTT[AAA]AC == AAA | add 1 to the kmer list

"""
# List comprehension solution
def count_kmer_listcomp(sequence, kmer):
    kmer_list = [sequence[pos:pos + len(kmer)] for pos in range(len(sequence) - (len(kmer) - 1))]
    return kmer_list.count(kmer)

"""
seq = "AATTAAAAC"
kmer = "AAA"
kmer_list <= add [AAT]TAAAAC
kmer_list <= add A[ATT]AAAAC
kmer_list <= add AA[TTA]AAAAC
...
count "AAA" in kmer_list
"""

# RegExp solution
def count_kmer_regexp(sequence, kmer):
    return len(findall(f'(?=({kmer}))',sequence))

# Test area

seq = "AAAGAAAAATTGA" * 100000
kmer = "AA"

start_time = perf_counter()
print("Loop k-mer count: ", end=' ')
print(count_kmer_loop(seq,kmer))
elapsed_time = perf_counter()
execution_time = elapsed_time - start_time
print(f'Loop took: {(elapsed_time - start_time):0.7f}s\n')

start_time = perf_counter()
print("List comprehension k-mer count: ", end=' ')
print(count_kmer_listcomp(seq,kmer))
elapsed_time = perf_counter()
execution_time = elapsed_time - start_time
print(f'List took: {(elapsed_time - start_time):0.7f}s\n')

start_time = perf_counter()
print("Regular expression k-mer count: ", end=' ')
print(count_kmer_regexp(seq,kmer))
elapsed_time = perf_counter()
execution_time = elapsed_time - start_time
print(f'Regexp took: {(elapsed_time - start_time):0.7f}s\n')

#start_time = perf_counter()

# call your code here

#elapsed_time = perf_counter()
#execution_time = elapsed_time - start_time