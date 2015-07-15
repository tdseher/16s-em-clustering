#!/bin/bash

# 16S clustering and assignment pipeline
# Copyright 2014 Thaddeus Seher

# If you use this 16S clustering and assignment pipeline, please
# cite our paper. Koch et al 2015

# Positional parameters:
#  1..n = FASTQ file for each sample

# Description
#  Creates a bunch of files and folders in the current working directory.
#  Aligns sequences to 16S database
#  Clusters sequences at specified threshold
#  Estimates abundance of each cluster in each sample
#  Assigns taxonomy to each cluster
#  Outputs a biom_table and biom file
#  Designed to be used with Qiime 1.8.0

# this old way doesn't work unless you call this program through its absolute path
#SCRIPT_PATH=$(dirname $0)
#SCRIPT_PATH=$(readlink -m $SCRIPT_PATH)

# this new way should work
SCRIPT_PATH=$(readlink -m $0)
SCRIPT_PATH=$(dirname $SCRIPT_PATH)
#SCRIPT_PATH=/data/control/programs/thad-16s-scripts

# Specify number of threads to use
THREADS=60
MAX_HOMOP=10
MIN_LENGTH=200
MAX_LENGTH=500
#ALIGNMENT_START=13862 # optimize
#ALIGNMENT_END=23444 # optimize
ALIGNMENT_START=optimize
ALIGNMENT_END=optimize
REFERENCE_ALIGNMENTS=/data/control/programs/mothur-1.33.3/referenceFiles/silva-release119-2.27-nr/silva.nr_v119.align
REFERENCE_NONCHIMERAS=/data/thaddeus/koch_projects/koch_miseq_192_samples/sel03n.masked.fixed.fasta
MIN_SIZE=10
IDENTITY=0.97
#DISTANCE=$(echo "100-(100*$IDENTITY)" | bc)
#DISTANCE=${DISTANCE%.*} # should be 3
DISTANCE=$(echo "1-$IDENTITY" | bc)
QIIME_ACTIVATE=/data/control/programs/qiime_software/activate.sh
QIIME_TAXONOMY=/data/control/programs/qiime_software/gg_otus-13_8-release/taxonomy/99_otu_taxonomy.txt
QIIME_REPSET=/data/control/programs/qiime_software/gg_otus-13_8-release/rep_set/99_otus.fasta
REVCOMP=1 # 1 or 0

## cat ../../sequences/filtered/Koch26S_Samples/filter/dmux/km.192.bcm.merged.Ch.2.* ../../sequences/filtered/Koch26S_Samples/filter/dmux/km.192.bcm.merged.CH.2.* > km.192.bcm.merged.ch2.fastq
#ALL_READS_FASTQ=$1
##ALL_READS_FASTA=$(echo $ALL_READS_FASTQ | perl -pe 's/\.fastq$/\.fasta/')
## Shift all positional parameters down by one
#shift;

# First, we concatenate the FASTQ files from each sample together into a single file
rm all.fastq
if [ "$REVCOMP" == 1 ]
then
  for SAMPLE_FASTQ in "$@"
  do
    python ${SCRIPT_PATH}/revcomp-fastq.py ${SAMPLE_FASTQ} >> all.fastq
  done
else
  cat "$@" > all.fastq
fi

# We use Mothur to convert to a FASTA file, filter sequences by length and max allowed homopolymer, and filter by alignment
# We don't use Mothur's pre.cluster() because it makes "usearch -cluster_otus" fail

if [ "$ALIGNMENT_START" == "optimize" ] && [ "$ALIGNMENT_END" == "optimize" ]
then
  mothur "#fastq.info(fastq=all.fastq, qfile=F); screen.seqs(maxhomop=${MAX_HOMOP}, minlength=${MIN_LENGTH}, maxlength=${MAX_LENGTH}); unique.seqs(); align.seqs(reference=${REFERENCE_ALIGNMENTS}, processors=${THREADS}); screen.seqs(optimize=start-end, criteria=90); filter.seqs(vertical=T, trump=.); deunique.seqs();" &>mothur.err
elif [ "$ALIGNMENT_START" == "optimize" ]
then
  mothur "#fastq.info(fastq=all.fastq, qfile=F); screen.seqs(maxhomop=${MAX_HOMOP}, minlength=${MIN_LENGTH}, maxlength=${MAX_LENGTH}); unique.seqs(); align.seqs(reference=${REFERENCE_ALIGNMENTS}, processors=${THREADS}); screen.seqs(optimize=start, criteria=90, end=${ALIGNMENT_END}); filter.seqs(vertical=T, trump=.); deunique.seqs();" &>mothur.err
elif [ "$ALIGNMENT_END" == "optimize" ]
then
  mothur "#fastq.info(fastq=all.fastq, qfile=F); screen.seqs(maxhomop=${MAX_HOMOP}, minlength=${MIN_LENGTH}, maxlength=${MAX_LENGTH}); unique.seqs(); align.seqs(reference=${REFERENCE_ALIGNMENTS}, processors=${THREADS}); screen.seqs(optimize=end, criteria=90, start=${ALIGNMENT_START}); filter.seqs(vertical=T, trump=.); deunique.seqs();" &>mothur.err
else
  mothur "#fastq.info(fastq=all.fastq, qfile=F); screen.seqs(maxhomop=${MAX_HOMOP}, minlength=${MIN_LENGTH}, maxlength=${MAX_LENGTH}); unique.seqs(); align.seqs(reference=${REFERENCE_ALIGNMENTS}, processors=${THREADS}); screen.seqs(start=${ALIGNMENT_START}, end=${ALIGNMENT_END}); filter.seqs(vertical=T, trump=.); deunique.seqs();" &>mothur.err
fi

# We get the filename of the FASTA file Mothur filtered for us
MOTHUR_OUTPUT_FASTA=$(tail -n 4 mothur.err | head -n 1)
MOTHUR_OUTPUT_PREFIX=$(echo $MOTHUR_OUTPUT_FASTA | perl -pe 's/\.fasta$//')

# We remove the gaps in the sequences
perl ${SCRIPT_PATH}/fasta_remove_gaps.pl ${MOTHUR_OUTPUT_PREFIX}.fasta > ${MOTHUR_OUTPUT_PREFIX}.nogaps.fasta

# We append a "size=1;" tag to the FASTA headers so we can feed the sequences into USearch
perl -e 'while (<>){chomp; my $line = $_; if (substr($line, 0, 1) eq ">") { print "$line;size=1;\n";}else{print "$line\n";}}' ${MOTHUR_OUTPUT_PREFIX}.nogaps.fasta > ${MOTHUR_OUTPUT_PREFIX}.nogaps.sized.fasta

# We dereplicate the sequences, appending the size=N tag to the FASTA header
#usearch -derep_fulllength ${MOTHUR_OUTPUT_PREFIX}.nogaps.fasta -sizeout -output ${MOTHUR_OUTPUT_PREFIX}.nogaps.derep.fasta
# We use a custom script because we don't own the 64-bit version of usearch
python ${SCRIPT_PATH}/dereplicate-v2.py ${MOTHUR_OUTPUT_PREFIX}.nogaps.fasta > ${MOTHUR_OUTPUT_PREFIX}.nogaps.derep.fasta

# We remove any sequences with size=1
usearch -sortbysize ${MOTHUR_OUTPUT_PREFIX}.nogaps.derep.fasta -output ${MOTHUR_OUTPUT_PREFIX}.nogaps.derep.min2.fasta -minsize 2

# We make our OTUs using USearch. This also performs an initial round of chimera checking.
#usearch -cluster_otus ${MOTHUR_OUTPUT_PREFIX}.sized.nohyphens.fasta -otu_radius_pct ${DISTANCE} -otus otus.fasta &> otus.err
usearch -cluster_otus ${MOTHUR_OUTPUT_PREFIX}.nogaps.derep.min2.fasta -id ${DISTANCE} -otus otus.fasta &> otus.err

# We set some convenience variables
QUERY=${MOTHUR_OUTPUT_PREFIX}.nogaps.sized.fasta
SUBJECT=otus.fasta
ALIGNMENT=alignment1

# We map our filtered reads to our OTUs, reporting all redundant alignments
usearch -usearch_global ${QUERY} -db ${SUBJECT} -strand plus -id ${IDENTITY} -maxaccepts 0 -maxrejects 0 -uc ${ALIGNMENT}.uc -uc_allhits -dbmatched ${ALIGNMENT}.dbm -sizeout

# We annotate the otus.fasta file with only the counts of reads that map to one OTU (we count only unique alignments)
python ${SCRIPT_PATH}/unique-uc2fasta.py otus.fasta ${ALIGNMENT}.uc > otus.unique.fasta

# We filter the OTUs, keeping only ones with at least $MIN_SIZE reads that map exclusively
usearch -sortbysize otus.unique.fasta -output otus.unique.sorted.fasta -minsize ${MIN_SIZE}

# convenience variables
SUBJECT=otus.unique.sorted.fasta
ALIGNMENT=alignment2
EXPRESS=express1
ABUNDANCE=otus

# We need to get the abundances of our OTUs, so we do the following
# We map our filtered reads to our OTUs again, reporting all redundant alignments
usearch -usearch_global ${QUERY} -db ${SUBJECT} -strand plus -id ${IDENTITY} -maxaccepts 0 -maxrejects 0 -uc ${ALIGNMENT}.uc -uc_allhits -dbmatched ${ALIGNMENT}.dbm -sizeout

# USearch version 8's beta doesn't output redundant *.sam output yet.
# So we convert our *.uc alignments to *.sam format
# We recommend using the "samtools view -H" command, which just prints the header,
# rather than the "samtools view -h" command, which prints the header and the entire *.sam file
# because CIGAR strings for some reason aren't properly parsed
rm ${ALIGNMENT}.sam
touch ${ALIGNMENT}.sam
samtools view -H -T ${SUBJECT} ${ALIGNMENT}.sam > ${ALIGNMENT}.sam

# uc2sam.py works alright, but for some reason SAM thinks its CIGAR strings are improperly-formatted
python ${SCRIPT_PATH}/uc2sam.py all.fastq ${ALIGNMENT}.uc >> ${ALIGNMENT}.sam

# no need to shuffle
# We shuffle the SAM file because eXpress expects a random sequence order
# and we append it to the header
#shuf ${ALIGNMENT}.preshuf.sam >> ${ALIGNMENT}.sam

# We remove the pre-shuffled SAM file to save space
#rm ${ALIGNMENT}.preshuf.sam

printf "%.0f" 2743410360.320

READ_LENGTH_MAX=$(cut -f 10 ${ALIGNMENT}.sam | perl -e 'my $max = 0; while(<>) { if (length($_) > $max) {$max = length($_);}} print "$max";')
READ_LENGTH_MEAN=$(grep -v '^@' ${ALIGNMENT}.sam | awk '!_[$1]++' | cut -f 10 | perl -e 'my $sum = 0; my $num = 0; while(<>) { $sum += length($_); $num += 1; } print $sum/$num;')
READ_LENGTH_MEAN=$(printf "%.0f" $READ_LENGTH_MEAN)
READ_LENGTH_STD=$(grep -v '^@' ${ALIGNMENT}.sam | awk '!_[$1]++' | cut -f 10 | perl -e 'use Math::NumberCruncher; my @data; while(<>) { push(@data, length($_)); } print Math::NumberCruncher::StandardDeviation(\@data);')
READ_LENGTH_STD=$(printf "%.0f" $READ_LENGTH_STD)

# We run eXpress to estimate the true sequence counts
express ${SUBJECT} ${ALIGNMENT}.sam --max-read-len ${READ_LENGTH_MAX} --no-bias-correct --frag-len-mean ${READ_LENGTH_MEAN} --frag-len-stddev ${READ_LENGTH_STD} --f-stranded --logtostderr --output-dir ${EXPRESS} 2> ${EXPRESS}.err

# We annotate the "size=N;" tag of the OTU file with the eXpress counts, and store in a *.fasta file
python ${SCRIPT_PATH}/express2fasta.py ${SUBJECT} ${EXPRESS}/results.xprs > ${ABUNDANCE}.abundance.fasta

# Now we do a round of de novo chimera filtering
usearch -uchime_denovo ${ABUNDANCE}.abundance.fasta -uchimeout otus.abundance.dn-uchime.results -uchimealns otus.abundance.dn-uchime.aln -chimeras otus.abundance.dn-uchime.chimeras.fasta -nonchimeras otus.abundance.dn-uchime.nonchimeras.fasta

# We do a reference-based chimera filter
usearch -uchime_ref otus.abundance.dn-uchime.nonchimeras.fasta -db ${REFERENCE_NONCHIMERAS} -strand plus -uchimeout otus.abundance.dn-uchime.nonchimeras.FTref.results -uchimealns otus.abundance.dn-uchime.nonchimeras.FTref.aln -chimeras otus.abundance.dn-uchime.nonchimeras.FTref.chimeras.fasta -nonchimeras otus.abundance.dn-uchime.nonchimeras.FTref.nonchimeras.fasta

# We map all our reads again, convert to SAM, estimate abundance, and do a length filter in order to get our final list of OTUs
QUERY=${MOTHUR_OUTPUT_PREFIX}.nogaps.sized.fasta
SUBJECT=otus.abundance.dn-uchime.nonchimeras.FTref.nonchimeras.fasta
ALIGNMENT=alignment3
EXPRESS=express2
ABUNDANCE=chimera-free-otus
usearch -usearch_global ${QUERY} -db ${SUBJECT} -strand plus -id ${IDENTITY} -maxaccepts 0 -maxrejects 0 -uc ${ALIGNMENT}.uc -uc_allhits -dbmatched ${ALIGNMENT}.dbm -sizeout
rm ${ALIGNMENT}.sam
touch ${ALIGNMENT}.sam
samtools view -H -T ${SUBJECT} ${ALIGNMENT}.sam > ${ALIGNMENT}.sam
python ${SCRIPT_PATH}/uc2sam.py all.fastq ${ALIGNMENT}.uc >> ${ALIGNMENT}.sam

READ_LENGTH_MAX=$(cut -f 10 ${ALIGNMENT}.sam | perl -e 'my $max = 0; while(<>) { if (length($_) > $max) {$max = length($_);}} print "$max";')
READ_LENGTH_MEAN=$(grep -v '^@' ${ALIGNMENT}.sam | awk '!_[$1]++' | cut -f 10 | perl -e 'my $sum = 0; my $num = 0; while(<>) { $sum += length($_); $num += 1; } print $sum/$num;')
READ_LENGTH_MEAN=$(printf "%.0f" $READ_LENGTH_MEAN)
READ_LENGTH_STD=$(grep -v '^@' ${ALIGNMENT}.sam | awk '!_[$1]++' | cut -f 10 | perl -e 'use Math::NumberCruncher; my @data; while(<>) { push(@data, length($_)); } print Math::NumberCruncher::StandardDeviation(\@data);')
READ_LENGTH_STD=$(printf "%.0f" $READ_LENGTH_STD)

express ${SUBJECT} ${ALIGNMENT}.sam --max-read-len ${READ_LENGTH_MAX} --no-bias-correct --frag-len-mean ${READ_LENGTH_MEAN} --frag-len-stddev ${READ_LENGTH_STD} --f-stranded --logtostderr --output-dir ${EXPRESS} 2> ${EXPRESS}.err
python ${SCRIPT_PATH}/express2fasta.py ${SUBJECT} ${EXPRESS}/results.xprs > ${ABUNDANCE}.abundance.fasta
usearch -sortbysize ${ABUNDANCE}.abundance.fasta -output ${ABUNDANCE}.abundance.sorted.fasta -minsize ${MIN_SIZE}

# Now we obtain the OTU counts from each sample
for SAMPLE_FASTQ in "$@"
do
  SAMPLE_DIR=${SAMPLE_FASTQ##*/}
  SAMPLE_DIR=${SAMPLE_DIR%.fastq}
  mkdir ${SAMPLE_DIR}
  
  # convenience vars
  SUBJECT=chimera-free-otus.abundance.sorted.fasta
  QUERY=${SAMPLE_DIR}.fasta
  ALIGNMENT=alignment
  
  if [ "$REVCOMP" == 1 ]
  then
    python ${SCRIPT_PATH}/revcomp-fastq.py ${SAMPLE_FASTQ} > ${SAMPLE_DIR}/${SAMPLE_DIR}.rc.fastq
    perl ${SCRIPT_PATH}/fastq2fasta.pl ${SAMPLE_DIR}/${SAMPLE_DIR}.rc.fastq > ${SAMPLE_DIR}/${QUERY}
  else
    perl ${SCRIPT_PATH}/fastq2fasta.pl ${SAMPLE_FASTQ} > ${SAMPLE_DIR}/${QUERY}
  fi
  
  usearch -usearch_global ${SAMPLE_DIR}/${QUERY} -db ${SUBJECT} -strand plus -id ${IDENTITY} -maxaccepts 0 -maxrejects 0 -uc ${SAMPLE_DIR}/${ALIGNMENT}.uc -uc_allhits -dbmatched ${SAMPLE_DIR}/${ALIGNMENT}.dbm -sizeout
  touch ${SAMPLE_DIR}/${ALIGNMENT}.sam
  samtools view -H -T ${SUBJECT} ${SAMPLE_DIR}/${ALIGNMENT}.sam > ${SAMPLE_DIR}/${ALIGNMENT}.sam
  python ${SCRIPT_PATH}/uc2sam.py ${SAMPLE_FASTQ} ${SAMPLE_DIR}/${ALIGNMENT}.uc >> ${SAMPLE_DIR}/${ALIGNMENT}.sam
  
  READ_LENGTH_MAX=$(cut -f 10 ${SAMPLE_DIR}/${ALIGNMENT}.sam | perl -e 'my $max = 0; while(<>) { if (length($_) > $max) {$max = length($_);}} print "$max";')
  READ_LENGTH_MEAN=$(grep -v '^@' ${SAMPLE_DIR}/${ALIGNMENT}.sam | awk '!_[$1]++' | cut -f 10 | perl -e 'my $sum = 0; my $num = 0; while(<>) { $sum += length($_); $num += 1; } print $sum/$num;')
  READ_LENGTH_MEAN=$(printf "%.0f" $READ_LENGTH_MEAN)
  READ_LENGTH_STD=$(grep -v '^@' ${SAMPLE_DIR}/${ALIGNMENT}.sam | awk '!_[$1]++' | cut -f 10 | perl -e 'use Math::NumberCruncher; my @data; while(<>) { push(@data, length($_)); } print Math::NumberCruncher::StandardDeviation(\@data);')
  READ_LENGTH_STD=$(printf "%.0f" $READ_LENGTH_STD)
  
  express ${SUBJECT} ${SAMPLE_DIR}/${ALIGNMENT}.sam --max-read-len ${READ_LENGTH_MAX} --no-bias-correct --frag-len-mean ${READ_LENGTH_MEAN} --frag-len-stddev ${READ_LENGTH_STD} --f-stranded --logtostderr --output-dir ${SAMPLE_DIR} 2> ${SAMPLE_DIR}/express.err
  python ${SCRIPT_PATH}/express2fasta.py ${SUBJECT} ${SAMPLE_DIR}/results.xprs > ${SAMPLE_DIR}/${SAMPLE_DIR}.abundance.fasta
  
  ln -s ${SAMPLE_DIR}/${SAMPLE_DIR}.abundance.fasta ${SAMPLE_DIR}.toadd
done

# Load Qiime 1.8.0 into the environment
source ${QIIME_ACTIVATE}

# Make symlinks to Qiime's 97% OTU *.fasta and taxonomy to a user-owned directory to get around the PyCogent permission error
mkdir qiime_files
ln -s ${QIIME_TAXONOMY} qiime_files/
ln -s ${QIIME_REPSET} qiime_files/

QIIME_TAXONOMY_USER=qiime_files/$(basename $QIIME_TAXONOMY)
QIIME_REPSET_USER=qiime_files/$(basename $QIIME_REPSET)

# Run Qiime on OTUs to identify species
assign_taxonomy.py -i chimera-free-otus.abundance.sorted.fasta -t ${QIIME_TAXONOMY_USER} -r ${QIIME_REPSET_USER} -m blast

# Now that we have assigned taxonomies to our OTU sequences, we convert all our data to a *.biom_table.txt file
# This script needs to be modified to take an argument that replaces the extraneous portions of the filename with empty strings
# (i.e. remove the "km.192.bcm.merged." from the beginning, and remove  ".toadd" from the tail of the filenames)
python ${SCRIPT_PATH}/qiime_assignments+usearch_mappings_to_otus2biom_table.py blast_assigned_taxonomy/chimera-free-otus.abundance.sorted_tax_assignments.txt *.toadd > experiment.biom_table.txt

# Some taxa have multiple 16S genes, so we merge their count and species information
python ${SCRIPT_PATH}/qiime_collapse_redundant_species_mappings2biom_table.py experiment.biom_table.txt > experiment.collapsed.biom_table.txt

# We convert the *.biom_table.txt file to a *.biom file
biom convert -i experiment.collapsed.biom_table.txt -o experiment.collapsed.dense.biom --table-type="otu table" --matrix-type=dense

# Fix the *.biom file so it can be analyzed with Qiime
python ${SCRIPT_PATH}/fix_qiime_biom_taxonomy.py experiment.collapsed.dense.biom > experiment.collapsed.dense.taxfix.biom

# Beautify the *.biom file so it is human-readable
python ${SCRIPT_PATH}/beautify_biom.py experiment.collapsed.dense.taxfix.biom > experiment.collapsed.dense.taxfix.beautified.biom

# Finished!
