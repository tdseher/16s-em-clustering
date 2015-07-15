#!/usr/bin/env python

# script to look at a *.uc file created with usearch
# and output only alignments that are unique (don't map to more than one subject)

# outputs a FASTA file with the size=N; tag

# input:
#  *.fasta file with sequences
#  *.uc file with alignments
# output:
#  *.fasta file with size=N; tag

import sys
import re
from collections import Counter

def load_fasta(fasta):
    flo = open(fasta, 'r')
    contig = {}
    name = None
    for line in flo:
        line = line.rstrip()
        if line.startswith('>'):
            name = re.split(r'\s+', line[1:], 1)[0]
            contig[name] = ''
        else:
            contig[name] += line
    flo.close()
    
    return contig

def load_uc(filename):
    flo = open(filename, 'r')
    same_query_lines = []
    output = []
    for line in flo:
        sline = line.rstrip().split("\t")
        if (len(same_query_lines) > 0):
            if (sline[8] == same_query_lines[0][8]):
                same_query_lines.append(sline)
            else:
                if (len(same_query_lines) == 1):
                    output += same_query_lines
                same_query_lines = [sline]
        else:
            same_query_lines.append(sline)
    
    if (len(same_query_lines) == 1):
        output += same_query_lines
    same_query_lines = []
    flo.close()
    return output

def print_unique(contigs, uc_lines):
    #targets = list(set(map(lambda x: x[9], uc_lines)))
    targets = Counter(map(lambda x: x[9], uc_lines))
    for header in contigs:
        #count = 0
        #for subject in map(lambda x: x[9], uc_lines):
        #for subject in targets:
        #    if (subject == header):
        #        count += 1
        
        #print ">" + header + ";size=" + str(count)
        new_header = re.sub(r'\s*;\s*size\s*=\s*\d+;?', '', header)
        print ">" + new_header + ";size=" + str(targets[header]) + ';'
        
        print contigs[header]

def check_arguments():
    if ((len(sys.argv) <= 2) or ('-h' in sys.argv[1:]) or ('--help' in sys.argv[1:])):
        print "USAGE python unique-uc2fasta.py input.fasta input.uc > output.fasta"
        return False
    return True

def main():
    if check_arguments():
        contigs = load_fasta(sys.argv[1])
        uc_lines = load_uc(sys.argv[2])
        print_unique(contigs, uc_lines)

if __name__ == '__main__':
    main()
