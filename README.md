# learning_python
This repo contains the work from MCB 198

# Library
seqlib.py is a library of essential functions that can be found throughout the following programs. Currently, it contains:

* read_fasta
	* reads a fasta file and returns each sequence name and sequence
* read_fastq
	* reads a fastq file and returns a list of 3 items: sequence ID, sequence, and quality
* gc
	* gives gc content of a provided sequence (as a string)
*gc_skew
	* calculates gc_skew of a provided sequence (as a string)
* Entropy
	* calculates shannons entropy for a provided sequence (as a string)

# Programs (completed)

## comparefeature.py
this program compares two feature files and searches for overlap. It takes 2 files as input, which contain 3 columns of data: chromosome ID, feature start bp coordinate, and feaure end bp coordinate. The user can also provide an overlap threshold, which filters the results to only include matches that overlap by the threshold amount or greater.

## entropy.py
this program calculates and masks regions of low entropy within a fasta file. it takes a file path and window size, with options for a threshold entropy level, step size, and a modification option for lower case modifications. The default entropy is 1. If a region of window size x has entropy below the threshold, the nucleotides in that window will be replaced with either "N" or "n", depending on the mod option. Entropy is calculated using the Shannon's entropy formula:

* H (total)= p*log2(p) + H
	* where p is the overall frequency of a nucleotide
	and H is the entropy

## qualityfastq.py
this program takes a fastq file, and returns the average quality score per nucleotide position. Output looks like the following:
* "Position 1, average: X"
	* where position indicates which nucleotide Position
	* and average gives the average quality score for that position.

## gcwindow.py
this program displays the gc content for a window of size x, given a fasta file and window size. It has options for a window step size, and a skew argument. If gc skew is selected, the program will also provide gc skew calculations for each windo using the following formula:
* skew=((g-c)/(g+c))
	* where g is the G count in the window, and c is the C count in the window.

## genomestats.py
genomestats.py provides a basic overview of a fasta file. It returns a report containing:
* total size of the sequence (bp length)
* number of contigs
* shortest and longest contig
* average and median contig size
* N50
* GC content
* nucleotide counts

## kmer.py
this program counts the kmers in a provided sequence. It takes in a path to a fasta file, and has optional arguments for kmer size (--k), reverse complement reading (--r), and a display option to only show x amount of kmers (--limit). The default kmer size is k=3.

* note: this program will not register kmers that contain "N"

## needlemanwunsch.py
this program simulates the Needleman Wunsch global alignment algorithm. It takes 2 files as well as match, mismatch, and indel point values, and returns a optimal sequence alignment.

## smithwaterman.py
this program simulates the Smith Waterman local alignment algorithm. It takes 2 files as well as match, mismatch, and indel point values, and returns a optimal sequence alignment.
