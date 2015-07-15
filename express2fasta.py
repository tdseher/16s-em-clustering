#!/usr/bin/env python

import sys
import re

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

def load_xprs(filename):
    #  0 bundle_id
    #  1 target_id <-- contig name
    #  2 length
    #  3 eff_length
    #  4 tot_counts
    #  5 uniq_counts
    #  6 est_counts <-- estimated count we want
    #  7 eff_counts
    #  8 ambig_distr_alpha
    #  9 ambig_distr_beta
    # 10 fpkm
    # 11 fpkm_conf_low
    # 12 fpkm_conf_high
    # 13 solvable
    # 14 tpm
    
    output = {}
    
    flo = open(filename, 'r')
    flo.next()
    for line in flo:
        sline = line.rstrip().split("\t")
        output[sline[1]] = int(round(float(sline[6])))
        
    flo.close()
    
    return output

def print_counts(contigs, counts):
    # prints the FASTA sequences in order from most-abundant to least-abundant
    for header in sorted(counts, key=lambda x: counts[x], reverse=True):
        new_header = re.sub(r'\s*;\s*size\s*=\s*\d+;?', '', header)
        print ">" + new_header + ";size=" + str(counts[header]) + ';'
        print contigs[header]

def check_arguments():
    if ((len(sys.argv) <= 2) or ('-h' in sys.argv[1:]) or ('--help' in sys.argv[1:])):
        print "USAGE python express2fasta.py input.fasta results.xprs > output.fasta"
        return False
    return True

def main():
    if check_arguments():
        contigs = load_fasta(sys.argv[1])
        counts = load_xprs(sys.argv[2])
        print_counts(contigs, counts)

if __name__ == '__main__':
    main()
