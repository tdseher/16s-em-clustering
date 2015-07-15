#!/usr/bin/env python

import sys
import argparse
import string

def revcomp(dna):
	complements = string.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
	return dna.translate(complements)[::-1]

parser = argparse.ArgumentParser(
	description="This script prints the reverse complement of the input FASTQ file. \
		Copyright 2013 Thaddeus Seher.",
	formatter_class = argparse.ArgumentDefaultsHelpFormatter
)
# Add mandatory arguments
parser.add_argument("fastq", help="FASTQ input file to parse")
# Add optional arguments
# none

args = parser.parse_args()

count = 0
flo = open(args.fastq, 'r')
for line in flo:
	line = line.rstrip()
	count += 1
	if ((count == 1) or (count == 3)):
		print line
	elif (count == 2):
		print revcomp(line)
	else:
		print line[::-1]
		count = 0
flo.close()
