import seqlib
import argparse
import sys

parser=argparse.ArgumentParser(description="Calculates and masks regions of low entropy")
parser.add_argument("--file", required=True, type = str, metavar = '<path>')
parser.add_argument("--window", required=True, type=int, metavar='<int>')
parser.add_argument("--threshold", required=False, default=1, type=int)
parser.add_argument("--mod", required=False, type=str)
parser.add_argument("--step", required=False, type=int, default=1)
parser.add_argument("--verbose", action="store_true")
arg = parser.parse_args()

def status(s):
    if arg.verbose:
        sys.stderr.write(s)
        sys.stderr.write("\n")

myfasta=seqlib.read_fasta(arg.file)
window=arg.window


for name, seq in myfasta:
    seq=seq.upper()
    i=0
    while (i < (len(seq)-window+1)):
        selection=seq[i:i+window]
        e=seqlib.entropy(selection)

        if e<arg.threshold:
            status("entering modifications")
            status(f"window:{selection} {e}")
            status(f"current seq {seq}")
            if arg.mod=="lower": selection=selection.lower()
            else:
                selection="N"*len(selection)
            newseq=seq.replace(seq[i:i+window], selection, 1)
            status(f"replacement seq: {newseq}")
            seq=newseq
            status(f"modified seq: {seq}")
            i+=arg.window
        else: i=i+arg.step
    print(">",end="")
    print(name)
    print(seq)
