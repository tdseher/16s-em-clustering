#!/usr/bin/env python

# reads in a qiime biom_table.txt
# writes to STDOUT a biom_table.txt

# if the 0th column is redundant, then add them together
import sys
flo = open(sys.argv[1], 'r')
otus = {}
for line in flo:
    line = line.rstrip()
    if line.startswith("#"):
        print line
    else:
        sline = line.split("\t")
        if (sline[0] in otus):
            for i in range(len(sline)):
                if (0 < i < len(sline) - 1):
                    otus[sline[0]][i] = str(int(otus[sline[0]][i]) + int(sline[i]))
        else:
            otus[sline[0]] = sline
flo.close()

for o in otus:
    print "\t".join(otus[o])