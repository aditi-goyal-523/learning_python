import seqlib
import sys


def fastq_quality(filename):
	myfile=seqlib.read_fastq(filename)
	qscores=[]
	for name, seq, quality in myfile:
		holder=[]
		for i in range(len(quality)):
			qs=quality[i]
			qv=ord(qs)-33
			p_qv=10**-(qv/10)
			holder.append(p_qv)
		qscores.append(holder)

	rows=len(qscores[0])
	avgs=[]
	for j in range(rows):
		col_total=0

		for i in range(len(qscores)):
			col_total+=qscores[i][j]

		avgs.append(col_total/(len(qscores)))

	return avgs

counter=0
for avgs in fastq_quality(sys.argv[1]):
    print("position", counter, "average:", avgs)
    counter+=1
