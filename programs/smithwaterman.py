import sys
import argparse
import seqlib

parser=argparse.ArgumentParser(description="Global alignment using SW alg")
parser.add_argument("--s1", required=True, type=str, metavar='<path>')
parser.add_argument("--s2", required=True, type=str, metavar='<path>')
parser.add_argument("--match", required=True, type=int)
parser.add_argument("--mismatch", required=True, type=int)
parser.add_argument("--indel", required=True, type=int)
arg = parser.parse_args()

'''
for debugging, this example is from wikipedia
seqa = "TGTTACGG"
seqb = "GGTTGACTA"
'''

file1= seqlib.read_fasta(arg.s1)
file2= seqlib.read_fasta(arg.s2)

seqa=""
seqb=""
for name, seq, in file1:
	seqa=seqa+seq

for name, seq in file2:
	seqb=seqb+seq
#make score and direction matrix

seqa=" "+seqa
seqb=" "+seqb
cols=len(seqa)
rows=len(seqb)
matrix = []
directions = []

	#setting up matrix dimensions
for i in range(rows):
	row= [None] * (cols)
	row_d=[None] * cols
	#set first number in each row = 0
	row[0]=0
	row_d[0]="N"
	#row.append(rest)
	matrix.append(row)
	directions.append(row_d)

#set first row to 0
matrix[0]= [0] * cols
directions[0]=["N"] * cols
#create a matrix printer
def printer(m):
	print("   ", "  ".join(seqa))
	for i in range(rows):
		print(seqb[i], "", m[i])

#scoring
match = arg.match
mismatch = arg.mismatch
indel = arg.indel

#make the direction key
dir_score={"D": 0, "T": None, "L": None, "N": 0}

'''
equations:
i=row
j=col
left=x[j-1]+indel
diag=x[i-1][j-1]+match/mismatch
top=x[i-1]+indel
chose the maximum of these 3 equations
'''

max_score=0
max_coord=(0, 0)

#begin scoring
for i in range(1, rows): #indexing the row
	for j in range(1, cols): #indexing the column
		if matrix[i][j] != None:
			continue
		else:
			dir_score["T"]=matrix[i - 1][j] + indel
			dir_score["L"]=matrix[i][j - 1] + indel
			#diagonal calculations
			if seqa[j] == seqb[i]:
				dir_score["D"]=matrix[i - 1][j - 1] + match
			else:
				dir_score["D"]=matrix[i - 1][j - 1] + mismatch

			#filling in the score matrix
			all_vals=dir_score.values()
			score=max(all_vals)
			matrix[i][j] = score

			#filling in direction matrix
			directions[i][j]=max(dir_score, key=dir_score.get)

			if score >= max_score:
				max_score = score
				max_coord = (i, j)

#traceback
#start at the max coordinate

i=i = max_coord[0]
j = max_coord[1]

seq_a=["SEQA"]
seq_b=["SEQB"]

while (i >= 0) and (j>=0):
	dir=directions[i][j]
	if dir == "D":
		seq_a.append(seqa[j])
		seq_b.append(seqb[i])
		i-=1
		j-=1
	elif dir == "L":
		seq_a.append(seqa[j])
		seq_b.append("-")
		j-=1
	elif dir == "T":
		seq_a.append("-")
		seq_b.append(seqb[i])
		i-=1
	elif dir == "N":
		print(seq_a)
		print(seq_b)
		break
printer(matrix)
