
import argparse
import sys

import seqlib

parser = argparse.ArgumentParser(
	description='Count k-mers in a sequence')
parser.add_argument('--fasta', required=True, type=str, metavar='<path>',
help='path to a fasta file, may be compressed')
parser.add_argument('--k', required=False, type=int, default=3,
 metavar="<int>", help='length of k-mer [%(default)i]')
parser.add_argument('--r', required=False, type=bool, default=False,
help='indicate if you want the reverse complement kmers by typing True')
parser.add_argument('--limit', required=False, type=int, metavar='<int>',
help='only show this many k-mers')
arg = parser.parse_args()


#counting
seq_count=0
nt_count=0
k_count = {}

#making the dictionary
for name, seq in seqlib.read_fasta(arg.fasta):
	seq = seq.upper()
	seq_count +=1
	nt_count += len(seq)
	for i in range(len(seq) - arg.k +1):
		kmer = seq[i:i+arg.k]
		if "N" not in kmer:
			if kmer not in k_count: k_count[kmer]=0
			k_count[kmer]+=1

#sorting by kmer count
sorted_kmers=sorted(k_count.items(), key=lambda x: x[1], reverse=True)
sorted_low=sorted(k_count.items(), key=lambda x: x[1])

#if limit is called, print only n of sorted counts
if arg.limit:
	n=0
	print("First", arg.limit, "sorted kmers:")
	for kmer, count in sorted_kmers:
		if n==arg.limit: break
		print(kmer, count)
		n+=1
	n=0
	print("\nLast", arg.limit, "sorted kmers:")
	for kmer, count in sorted_low:
		if n==arg.limit: break
		print(kmer, count)
		n+=1

else: print("sorted kmers:\n", sorted_kmers, "\n")


#making the reverse dictionary
if arg.r:
	keys=list(k_count.keys())
	counts=list(k_count.values())
	reverse=[]
	for str in keys:
		reverse.append(str[::-1])
		reverse_count=dict(zip(reverse, counts))
		reverse_sorted=sorted(reverse_count.items(), key = lambda x: x[1],
		reverse=True)
		reverse_sorted_low=sorted(reverse_count.items(), key = lambda x: x[1])

	#printing the reverse dictionary
	if arg.limit:
		n=0
		print("First", arg.limit, "sorted reverse complement kmers:")
		for kmer, count in reverse_sorted:
			if n==arg.limit: break
			print(kmer, count)
			n+=1
		n=0
		print("\nLast", arg.limit, "sorted reverse complement kmers:")
		for kmer, count in reverse_sorted_low:
			if n==arg.limit: break
			print(kmer, count)
			n+=1
	else: print("sorted reverse complement kmers:\n", sorted_kmers, "\n")
