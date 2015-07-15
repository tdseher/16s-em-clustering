#!/usr/bin/env python

import sys

counts = {}

flo = open(sys.argv[1], 'r')

header = None
sequence = ''

for line in flo:
    line = line.rstrip()
    if line.startswith(">"):
        if header:
            try:
                counts[sequence][1] += 1
            except KeyError:
                counts[sequence] = [header, 1]
        header = line
        sequence = ''
    else:
        sequence += line

if header:
    try:
        counts[sequence][1] += 1
    except KeyError:
        counts[sequence] = [header, 1]

for s in sorted(counts, key=lambda x: counts[x][1], reverse=True):
    if (counts[s][0][-1] != ';'):
        counts[s][0] += ';'
    print counts[s][0] + 'size=' + str(counts[s][1]) + ';'
    print s