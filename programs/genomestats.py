import seqlib
import argparse
import sys

parser = argparse.ArgumentParser(description="Displays stats about a genome")
parser.add_argument("--file", required=True, type=str, metavar='<path>')
arg = parser.parse_args()

myfasta = seqlib.read_fasta(arg.file)

contigs = 0
contig_lengths = []
nt_length = 0
nt_counts = {}

for name, seq in myfasta:
    seq = seq.upper()
    length = len(seq)
    contigs += 1
    contig_lengths.append(length)
    nt_length += length

    for i in range(length):
        n = seq[i]
        if n not in nt_counts: nt_counts[n] = 0
        nt_counts[n] += 1

avg = nt_length / contigs
contig_lengths.sort(reverse=True)

short = contig_lengths[-1]
long = contig_lengths[0]

if (contigs % 2 == 0): median = contig_lengths[(int((contigs / 2)) - 1)]
else: median = contig_lengths[int(((contigs - 1) / 2))]

n50_target = nt_length / 2

n50 = 0
counter = 0

for i in range(len(contig_lengths)):
    if counter < n50_target:
        counter += contig_lengths[i]
    else:
        n50 = contig_lengths[i - 1]
        break

c_count = nt_counts["C"]
g_count = nt_counts["G"]
a_count = nt_counts["A"]
t_count = nt_counts["T"]

gc = (c_count + g_count) / (c_count + g_count + a_count + t_count)

print("Final Report:")
print("Total size: ", nt_length)
print("Number of contigs: ", contigs)
print("Shortest contig: ", short)
print("Largest contig: ", long)
print("Average contig size: ", avg)
print("Median contig size: ", median)
print("N50: ", n50)
print("GC fraction: ", gc)
print("counts of each letter: ", nt_counts)
