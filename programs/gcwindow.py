import seqlib
import argparse
import sys

parser = argparse.ArgumentParser(description="Displays gc content for windows")
parser.add_argument("--file", required=True, type=str, metavar='<path>')
parser.add_argument("--size", required=True, type=int, metavar='<int>')
parser.add_argument("--step", required=False, type=int, default=1,
metavar='<int>')
parser.add_argument("--skew", required=False, type=bool)

arg=parser.parse_args()

myfasta=seqlib.read_fasta(arg.file)
window=arg.size

results=[]

for name, seq in myfasta:
	seq=seq.upper()
	i=0
	while (i < (len(seq)-window+1)):
		selection=seq[i:i+window]
		if arg.skew: results.append(list(name, i, seqlib.gc_skew(selection)))
		else:        results.append(list(name, i, seqlib.gc(selection)))
		i+=arg.step

print(results)
