#!/usr/bin/env python

import json
import biom
import sys


def dict_by_line(in_dict, indent):
	out_str = ""
	length = len(in_dict.keys())
	count = 0
	for k,v in in_dict.iteritems():
		count += 1
		if (count < length):
			out_str += '{0}{1}: {2},\n'.format(indent, json.dumps(k),json.dumps(v))
			#out_str += '"{0}": "{1}",\n'.format(k,v)
		else:
			out_str += '{0}{1}: {2},'.format(indent, json.dumps(k),json.dumps(v))
		
	return out_str
	

def list_by_line(in_list, indent):
	out_str = ""
	length = len(in_list)
	count = 0
	for v in in_list:
		count += 1
		if (count < length):
			out_str += '{0}{1},\n'.format(indent, json.dumps(v))
		else:
			out_str += '{0}{1}'.format(indent, json.dumps(v))
		
	return out_str


# open the .biom filename specified on the command line
filename = sys.argv[1]

# Read the .biom file into the data variable
with open(filename, 'r') as infile:
	all_data = json.load(infile)

# alternative:
#print json.dumps(data, sort_keys=True, indent=4)


rows = all_data.pop('rows')
columns = all_data.pop('columns')
data = all_data.pop('data')

print '{'
print dict_by_line(all_data, '  ')
print '  "rows": ['
print list_by_line(rows, '    ')
print '  ],'
print '  "columns": ['
print list_by_line(columns, '    ')
print '  ],'
print '  "data": ['
print list_by_line(data, '    ')
print '  ]'
print '}'


