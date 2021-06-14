import sys
import argparse
'''
parser=argparse.ArgumentParser(description="Global alignment using NW alg")
add arg for the file names (how to do with no dashes? find positional arg)
parser.add_argument("--match", required=True, type=int)
parser.add_argument("--mismatch", required=True, type=int)
parser.add_argument("--indel", required=True, type=nt)
arg = parser.parse_args()


'''

seqa = "GCATGGU"
seqb = "GATTACA"


def print_matrix(m):
    for i in m:
        print(i)
    print("---------------")
    

#1 constructing the grid
matrix = []
trace = []

for i in range(len(seqb)):
    row = []
    row_trace = []
    for j in range(len(seqa) - 1):
        row.append(None)
        row_trace.append(None)
    matrix.append(row)
    trace.append(row_trace)

#scoring
match = 1
#match=arg.match
mismatch = -1
#mismatch=arg.mismatch
indel = -1
#indel=arg.indel
'''
equations:
i=row
j=col
left=x[j-1]+indel
diag=x[i-1][j-1]+match/mismatch
top=x[i-1]+indel

chose the maximum of these 3 equations
'''
matrix[1][1] = 0
for j in range(2, len(matrix[1])):
    matrix[1][j] = (j - 1) * indel

for i in range(2, len(matrix)):
    matrix[i][1] = (i - 1) * indel

print_matrix(matrix)
for i in range(2, len(matrix)):
    for j in range(2, len(matrix[i])):
        #now they're all fair game
        if matrix[i][j] != None:
            continue
        else:
            print(matrix[i - 1][j])
            print(matrix[i][j - 1])
            top = (matrix[i - 1][j] + indel)
            left = (matrix[i][j - 1] + indel)

            #diagonal calculations
            row_nt = matrix[0][j]
            col_nt = matrix[i][0]
            if row_nt == col_nt: diag = (matrix[i - 1][j - 1] + match)
            else: diag = (matrix[i - 1][j - 1] + mismatch)
            if diag >= top and diag >= left:
                matrix[i][j] = diag
                trace[i][j] = "d"
            elif left >= top and left >= diag:
                matrix[i][j] = left
                trace[i][j] = "l"
            else:
                matrix[i][j] = top
                trace[i][j] = "t"

print_matrix(matrix)
print_matrix(trace)
#alignment
alignment = []
topseq = [matrix[0][len(matrix[0]) - 1]]
leftseq = [matrix[len(matrix) - 1][0]]

i = len(matrix) - 1
j = len(matrix[i]) - 1

while i > 1 and j > 1:
    pos = matrix[i][j]
    if trace[i][j] == "d":
        topseq.append(seqa[j - 1])
        leftseq.append(seqb[i - 1])
    elif trace[i][j] == "t":
        topseq.append("-")
        leftseq.append(seqb[i - 1])
    else:
        topseq.append(seqa[j - 1])
        leftseq.append("-")

print_matrix(matrix)
