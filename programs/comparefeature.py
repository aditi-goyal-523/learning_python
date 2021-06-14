import sys
import argparse

parser = argparse.ArgumentParser(description="compares features and searches for overlap")
parser.add_argument("--threshold", required=False, type = int, default=1, metavar = '<int>')
parser.add_argument("--f1", required=True, type = str, metavar = '<path>')
parser.add_argument("--f2", required=True, type = str, metavar = '<path>')
arg = parser.parse_args()

def read_features(file):
	features={}
	with open(file) as fp:
		for line in fp.readlines():
			(chr,b,e)=line.split()
			if chr not in features: features[chr]=[]
			features[chr].append((int(b), int(e)))
	return features

def overlap(a, b, t):
	b1, e1=a
	b2, e2=b

	if b1>=b2 and b1<=e2:
		if e2-b1>=t-1:
			return True

	elif e1>=b2 and e1<=e2:
		if e1-b2>=t-1:
			return True

	elif b1<=b2 and e1>=e2: return True
	elif b2<=b1 and e1<=e2: return True
	else: return False

f1=read_features(arg.f1)
f2=read_features(arg.f2)

overlaps=[]

for chr in f1:
	if chr not in f2:
		continue
	for f in f1[chr]:

		for g in f2[chr]:
			if overlap(f, g, arg.threshold):
				overlaps.append([chr, f, g])

print(len(overlaps))
