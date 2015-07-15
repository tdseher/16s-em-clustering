#!/usr/bin/env python

# fixes the taxonomy field for qiime *.biom files created from *.biom_table files
import json
import sys
import re


# open the .biom filename specified on the command line
filename = sys.argv[1]

# Read the .biom file into the data variable
with open(filename, 'r') as infile:
    all_data = json.load(infile)


# Determine the data type. Can be 'Taxonomy', 'Subsystems', or 'KO'
data_type = all_data['rows'][0]['metadata'].keys()[0]

for i in range(len(all_data['rows'])):
    line = all_data['rows'][i]['metadata'][data_type]
    sline = re.split(r'\s*;\s*', line)
    all_data['rows'][i]['metadata'][data_type] = sline
    #print "\t".join(line) + "\t" + "\t".join(map(str,data['data'][i]))

# rows = all_data.pop('rows')
# for r in rows:
#     r["metadata"]["taxonomy"]


print json.dumps(all_data)