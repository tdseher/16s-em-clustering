#!/usr/bin/env python

# Input:
#  the *.uc file from using usearch_global

# Field	Description
# 1	Record type S, H, C or N (see table below).
# 2	Cluster number (0-based).
# 3	Sequence length (S, N and H) or cluster size (C).
# 4	For H records, percent identity with target.
# 5	For H records, the strand: + or - for nucleotides, . for proteins.
# 6	Not used, parsers should ignore this field. Included for backwards compatibility.
# 7	Not used, parsers should ignore this field. Included for backwards compatibility.
# 8	Compressed alignment or the symbol '=' (equals sign). The = indicates that the query is 100% identical to the target sequence (field 10).
# 9	Label of query sequence (always present).
# 10	Label of target sequence (H records only).

# Record 	Description
# H	Hit. Represents a query-target alignment. For clustering, indicates the cluster assignment for the query. If maxaccepts > 1, only there is only one H record giving the best hit. To get the other accepts, use another type of output file, or use the -uc_allhits option (requires version 6.0.217 or later).
# S	Centroid (clustering only). There is one S record for each cluster, this gives the centroid (representative) sequence label in the 9th field. Redundant with the C record; provided for backwards compatibility.
# C	Cluster record (clustering only). The 3rd field is set to the cluster size (number of sequences in the cluster) and the 9th field is set to the label of the centroid sequence.
# N	No hit (for database search without clustering only). Indicates that no accepts were found. In the case of clustering, a query with no hits becomes the centroid of a new cluster and generates an S record instead of an N record.

# A compressed alignments represents an alignment in a compact format that does
# not include the sequence letters. The representation uses run-length encoding,
# as follows. Each column in the alignment is classified as M, D or I.
#
# Column	Description
# M	Match. A pair of letters.
# D	Delete. A gap in the target.
# I	Insert. A gap in the query.

# If there are n consecutive columns of type C, this is represented as nC.
# For example, 123M is 123 consecutive matches. As a special case, if n=1 then
# n is omitted. So for example, D5M2I3M represents an alignment of this form:
#   Query    XXXXXX--XXX
#   Target   -XXXXXXXXXX
#   Column   DMMMMMIIMMM


# Output:
#  a *.sam file

# Field	Description
# 1	Query sequence label (typically, a read label).
# 2	Integer representing the sum of integer flag values: 4=no hit, 16=rev-comp'd, 256=secondary hit (all hits except one, the top hit, have this flag).
# 3	Target sequence label (reference database label).
# 4	1-based position in target.
# 5	Mapping quality (ignored / set to * by USEARCH).
# 6	CIGAR string.
# 7	Ref name of mate (ignored / set to * by USEARCH).
# 8	Position of mate (ignored / set to * by USEARCH).
# 9	Target sequence length.
# 10	Query sequence. If soft clipping is used, this is the full-length query. If hard clipping is used, this is just the alignment segment of the query. The sam_softclip option specifies soft clipping. Default is hard clipping, because with soft clipping SAM files can be very large with long query sequences.
# 11	ASCII string with Phred scores (ignored / set to * by USEARCH).
# 12, 13..	Tag fields (optional). Currently, USEARCH generates the following tags (this may be subject to change): AS, XN, XM, XO, XG, NM, MD and YT. Most tags are ignored on input. MD is required in an input file in order to generate complete alignments if the database is not provided (see alnout option to sam_filter).

import sys
import re
def check_arguments():
    if ((len(sys.argv) <= 2) or ('-h' in sys.argv[1:]) or ('--help' in sys.argv[1:])):
        print "USAGE:"
        print "First, run the following:"
        print "usearch -usearch_global input.fasta -db otus.fasta -strand plus -id 0.97 -maxaccepts 0 -maxrejects 0 -uc input.uc -uc_allhits"
        print "Then run this:"
        print "python uc2sam.py input.fastq input.uc > output.sam"
        return False
    return True

def uc2cigar(cigar, regex):
    number_added = 0
    for m in regex.finditer(cigar):
        if len(m.group()) > 1:
            position = m.start() + number_added + len(m.group()) - 1
        else:
            position = m.start() + number_added + len(m.group())
        cigar = cigar[:position] + '1' + cigar[position:]
        number_added += 1
    return cigar

def uc2sam(uc_sline, regex):
    sam_sline = [None, '4', '*', '0', '0', '*', '*', '0', '0', None, None]
    
    sam_sline[0] = uc_sline[8] # query ID
    if (uc_sline[9] != '*'):
        sam_sline[1] = '0' # FLAG
        sam_sline[2] = uc_sline[9] # subject ID
        sam_sline[3] = '1' # subject mapping position
        sam_sline[4] = '255' # mapping quality
        
        #sam_sline[5] = uc_sline[7] # CIGAR string
        sam_sline[5] = uc2cigar(uc_sline[7], regex) # CIGAR string
        
        #sam_sline[6] = '*' # mate subject ID
        #sam_sline[7] = '0' # mate mapping position
        #sam_sline[8] = '0' # insert size
    
    return sam_sline

def uc_lines(uc_slines, regex):
    # 0: 'read mapped',
    # 4: 'read unmapped',
    # 256: 'not primary alignment',
    # 2048: 'supplementary alignment',
    
    # XT:A:U # tag for unique alignment
    # XT:A:R # tag for non-unique alignment
    # NM:i:0 # tag for aligned read edit distance
    # MD:Z:107G118T26 # tag for mismatch positions/bases
    # X0:i:1 # tag with number of best hits
    # X1:i:0 # number of non-best hits (suboptimal hits or secondary alignments)
    # XN:i:0 # tag for number of ambiguous bases (N's) if the XN tag is greater than zero
    # XO:i:0 # number of gap opens
    # XG:i:0 # number of gap extensions
    # XA # Alternative hits; format: (chr,pos,CIGAR,NM;)*
    # XS = suboptimal alignment score
    
    sam_lines = []
    
    for uc_sline in uc_slines:
        sam_sline = uc2sam(uc_sline, regex)
        if(len(uc_slines) == 1):
            #sam_sline[1] = '0' # FLAG
            #sam_sline[4] = '255' # mapping quality
            sam_sline.append('XT:A:U') # tag for unique alignment
        else:
            #sam_sline[1] = '0' # FLAG
            #sam_sline[4] = '255' # mapping quality
            sam_sline.append('XT:A:R') # tag for unique alignment
        
        if (sam_sline[1] == '4'):
            sam_sline.append('XM:i:0') # tag for unmapped reads, indicating the number of valid alignments suppressed
            # or XM:i:0	= Number of mismatches in the alignment?
        else:
            sam_sline.append('XA:i:0') # tag for the number of mismatches in the "seed" region of the alignment
        
        sam_lines.append(sam_sline)
    
    return sam_lines

def load_uc(filename):
    #regex = re.compile(r'^[D]|[MID][D]')
    regex = re.compile(r'^[ID]|[MID][ID]')
    slines = []
    same_query_lines = []
    flo = open(filename, 'r')
    for line in flo:
        uc_sline = line.rstrip().split("\t")
        if (len(same_query_lines) > 0):
            if (uc_sline[8] == same_query_lines[0][8]):
                same_query_lines.append(uc_sline)
            else:
                slines += uc_lines(same_query_lines, regex)
                same_query_lines = [uc_sline]
        else:
            same_query_lines.append(uc_sline)
    
    slines += uc_lines(same_query_lines, regex)
    same_query_lines = []
    
    flo.close()
    return slines

def load_fastq(filename):
    seqs = {}
    quals = {}
    
    rec = re.compile(r'^@([\S]*).*$')
    
    quad = []
    flo = open(filename, 'r')
    for line in flo:
        quad.append(line.rstrip())
        if (len(quad) == 4):
            m = rec.match(quad[0])
            header = m.group(1).replace(':', '_')
            seqs[header] = quad[1]
            quals[header] = quad[3]
            quad = []
    flo.close()
    
    return seqs, quals

def main():
    if check_arguments():
        seqs, quals = load_fastq(sys.argv[1])
        slines = load_uc(sys.argv[2])
        for i in range(len(slines)):
            new_header = re.split(r'\s+', slines[i][0], 1)[0]
            new_header = re.sub(r'\s*;\s*size\s*=\s*\d+;?', '', new_header)
            new_header = new_header.replace(':', '_')
            slines[i][9] = seqs[new_header]
            slines[i][10] = quals[new_header]
            #print "\t".join(map(str, slines[i]))
            print "\t".join(slines[i])

if __name__ == '__main__':
    main()