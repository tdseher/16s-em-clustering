#!/usr/bin/env python

# load in the qiime classification file (blast_assigned_taxonomy/usearch.otus.dbm.sorted.dn-uchime.nonchimeras.FTref.nonchimeras.usearch.dbm.sorted_tax_assignments.txt)
# make a dict
# iterate through the lines
#  parse the 0th column to get the header
#  store taxonomy (col 1) and qiime id (col 3) in the dict

# load in a list of usearch-mapped FASTA files

# for each FASTA file
#  store the file name
#  sample_name = parsed file name
#  new dict
#  iterate through lines in file
#   store the FASTA header + size tag in dict

# print header row: #<OTU ID> <sample1> <sample2> ... <sampleN> <taxonomy>
# print each OTU on a separate row, and the counts in the columns

import sys
import re

def load_fasta(fasta):
    flo = open(fasta, 'r')
    contig = {}
    name = None
    for line in flo:
        line = line.rstrip()
        if line.startswith('>'):
            name = line[1:]
            contig[name] = ''
        else:
            contig[name] += line
    flo.close()
    
    return contig

# iterate through the qiime-blast matches
taxonomies = {}
otus = {}
assignment_flo = open(sys.argv[1], 'r')
for line in assignment_flo:
    m = re.match(r'^\s*#.*', line)
    if not m:
        sline = line.rstrip().split("\t")
        header = sline[0].split(";")[0] # strip the ';size=N;' tag
        taxonomies[header] = sline[1]
        otus[header] = sline[3]
assignment_flo.close()


# iterate through the usearch-mapped FASTA files
samples = {}
all_contigs = {}

for infile in sys.argv[2:]: # open the arguments as FLOs
    contigs = load_fasta(infile)
    
    name = infile.replace('km.192.bcm.merged.', '').replace('.toadd', '')
    samples[name] = {}
    
    for full_header in contigs:
        split_header = full_header.split(';')
        #header = split_header[0]
        #size = split_header[1]
        
        # update the list of all contigs
        all_contigs[split_header[0]] = contigs[full_header]
        
        samples[name][split_header[0]] = split_header[1].split('=')[1]

# print the table
print '# OTU ID\t' + '\t'.join(sorted(samples)) + '\ttaxonomy'
for otu in all_contigs:
    counts = []
    for s in sorted(samples):
        try:
            counts.append(samples[s][otu])
        except KeyError:
            counts.append('0')
    print "\t".join([otus[otu]] + counts + [taxonomies[otu]])
    
    