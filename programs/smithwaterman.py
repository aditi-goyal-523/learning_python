seqa = "TGTTACGG"
seqb = "GGTTGACTA"

#wewrite with traceback


def print_matrix(m):
    for i in m:
        print(i)
    print("---------------")


#1 constructing the grid
matrix = []
directions = []
header = list(seqa)
header.insert(0, None)
header.insert(0, None)
matrix.append(header)

blank = [None] * len(header)
matrix.append(blank)

for i in range(len(seqb)):
    row = []
    row.append(seqb[i])
    for j in range(len(header) - 1):
        row.append(None)
    matrix.append(row)

#scoring
match = 3
mismatch = -3
indel = -2
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

#start with matrix row 1 (i=1)
#start with matrix col 2 (j=2)
max_score = 0
max_coord = (1, 1)
for i in range(1, len(matrix)):
    for j in range(1, len(matrix[i])):
        if i == 1:
            #top and diagonal can't be calculated, just assign the left value
            if matrix[i][j] != None:
                continue
            else:
                matrix[i][j] = 0
        elif j == 1:
            #diagnonal and left can't be calculated, just assign top value
            if matrix[i][j] != None:
                continue
            else:
                matrix[i][j] = 0
        else:
            #now they're all fair game
            if matrix[i][j] != None:
                continue
            else:
                top = (matrix[i - 1][j] + indel)
                left = (matrix[i][j - 1] + indel)

                #diagonal calculations
                row_nt = matrix[0][j]
                col_nt = matrix[i][0]
                if row_nt == col_nt:
                    diag = (matrix[i - 1][j - 1] + match)
                else:
                    diag = (matrix[i - 1][j - 1] + mismatch)

                maximum = max(0, top, left, diag)
                matrix[i][j] = maximum
                if maximum >= max_score:
                    max_score = maximum
                    max_coord = (i, j)

#finding local maximum
maxvalue = matrix[1][1]
maxcoord = (1, 1)

for i in range(1, len(matrix)):
    for j in range(1, len(matrix[i])):
        valholder = matrix[i][j]
        if valholder > maxvalue:
            maxvalue = matrix[i][j]
            maxcoord = (i, j)
            print(f"new maximum set at {maxvalue} found at {maxcoord}")

#alignment
alignment = []
topseq = [matrix[0][maxcoord[1]]]
leftseq = [matrix[maxcoord[0]][0]]

i = maxcoord[0]
j = maxcoord[1]
pos = matrix[i][j]

while i > 1:
    '''
	for each position we basically need to redo the calculations from before
	use that to find directions
	then backtrace
	'''
    pos = matrix[i][j]
    top = (matrix[i - 1][j] + indel)
    left = (matrix[i][j - 1] + indel)

    #diagonal calculations
    row_nt = matrix[0][j]
    col_nt = matrix[i][0]
    if row_nt == col_nt: diag = (matrix[i - 1][j - 1] + match)
    else: diag = (matrix[i - 1][j - 1] + mismatch)
    if pos == diag:
        if matrix[i - 1][j - 1] != 0:
            print(pos)
            leftseq.insert(0, matrix[i - 1][0])
            topseq.insert(0, matrix[0][j - 1])
            print("T:", topseq)
            print("L:", leftseq)
            print("-------")
            i -= 1
            j -= 1
    elif pos == top:
        if matrix[i - 1][j] != 0:
            print(pos)
            leftseq.insert(0, matrix[i - 1][0])
            topseq.insert(0, "-")
            print("T:", topseq)
            print("L:", leftseq)
            print("-------")
            i -= 1
    elif pos == left:
        if matrix[i][j - 1] != 0:
            print(pos)
            leftseq.insert(0, "-")
            topseq.insert(0, matrix[0][j - 1])
            print("T:", topseq)
            print("L:", leftseq)
            print("-------")
            j -= 1
            continue
