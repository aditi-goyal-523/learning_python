import argparse
import sys

import seqlib

parser = argparse.ArgumentParser(
	description='Count k-mers in a sequence')
parser.add_argument('--fasta', required=True, type=str, metavar='<path>', help='path to a fasta file, may be compressed')
parser.add_argument('--k', required=False, type=int, default=3, metavar="<int.", help='length of k-mer [%(default)i]')
parser.add_argument('--limit', required=False, type=int, metavar='<int>', help='only show this many k-mers')
parser.add_argument('--verbose', action='store_true', help='print some diagnostic messages to stderr')
arge = parser.parse_args()

if arg.verbose:
	sys.stderr.write(f'Reading {arg.fasta}\n')
	
#counting
seq_count=0
nt_count=0
k_count= {}

for name, seq in seqlib.read_fasta(arg.fasta):
	seq_count +=1
	nt_count += len(seq)
	for i in range(len(seq) - arg.k +1):
		kmer = seq[i:i+arg.k]
		ifkmer not in k_count: k_count[kmer]=0
		k_count[kmer]+= 1
		
if arg.verbose:
	sys.stderr.write(f'{seq-count} sequences\n')
	sys.stderr.write(f'{nt_count} letters\n')
	sys.stderr.write(f'{len(k_count)} kmers\n')
	
#output
n=0
for kmer, count in sorted(k_count.items(), key=lambda itme: item[1], reverse = True):
	if arg.limit and n==arg.limit: break
	print(kmer, count)
	n+=1

