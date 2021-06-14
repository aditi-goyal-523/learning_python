import gzip
import math
import sys

def read_fasta(filename):
	name = None
	seqs = []

	fp = None
	if filename.endswith('.gz'): fp =gzip.open(filename, 'rt')
	elif filename=="-": fp=sys.stdin
	else: 	fp = open(filename)
	#---------------------------------------------------------------------
	for line in fp.readlines(): #for line in set of lines
		line = line.rstrip() 		#format the line so it is readable
		#-------------------------------------------------------------------
		if line.startswith('>'): 	#if the line starts with the > sign
			if len(seqs)>0:				#if there are already sequences in the list
				seq = "".join(seqs)			#seq is all the sequences already defined seperated by a space (i 												think this is if there are like...fragments under the same name? 												So it concatenates it all into one line)
				yield(name, seq)			#return the name, and the sequence
				name = line[1:]				#name holder is the line, first element of the tuple
				seqs = []					#now seqs is empty again
			else:						#else if the sequences is empty
				name = line[1:]				#the name is just the first element of that tuple. This is standard
		else:						#else if the line does not start with the > (meaning it is the sequence)
			seqs.append(line)			#just join that line to the seqence
	yield(name, ''.join(seqs))	#return the names, and the concatenated sequences
	fp.close()					#all don

def read_fastq(filename):
	name=None
	seqs=[]
	quality=None

	fp=None
	if filename.endswith('.gz'): fp=gzip.open(filename, 'rt')
	else: fp=open(filename)
	lines=[]
	for line in fp.readlines():
		lines.append(line)
	fp.close()
	for i in range(0, len(lines)-3, 4):
		name=lines[i][1:]
		seq=lines[i+1]
		quality=lines[i+3]
		name=name.rstrip()
		seq=seq.rstrip()
		quality=quality.rstrip()
		yield name, seq, quality

def gc(mysequence):
	mysequence=mysequence.upper()
	g=mysequence.count("G")
	c=mysequence.count("C")
	a=mysequence.count("A")
	t=mysequence.count("T")
	total = g+c+a+t
	if total == 0: return None
	gc=(g+c)/(total)
	return gc

def gc_skew(mysequence):
	mysequence=mysequence.upper()
	g=mysequence.count("G")
	c=mysequence.count("C")
	if (g+c) == 0: return None
	skew=((g-c)/(g+c))
	return skew

def entropy(seq):
	a=seq.count("A")
	c=seq.count("C")
	g=seq.count("G")
	t=seq.count("T")
	total=a+g+c+t
	p=(a/total, c/total, g/total, t/total)

	sum=0
	for pi in p: sum+=pi
	assert(math.isclose(sum, 1))

	H=0
	for pi in p:
		if pi==0: continue
		H+=pi*math.log2(pi)

	return -H
