#!/bin/sh

echo "***PAIRED END UNIFIED CONSENSUS STATISTICS***"
echo " "

echo "-------------RAW READS-------------"
samtools flagstat $1.temp.sort.bam
echo " "

echo "-------------SSCS-------------"
samtools flagstat $1_mem.sscs.sort.bam
echo " "

echo "-------------DCS-------------"
samtools flagstat $1_mem.dcs.sort.bam
echo " "

echo "-------------DCS FILTERING-------------"

echo "$(samtools view $1.dcs.filt.clipped.30end.realign.bam 17:7572523-7580325 | wc -l) Reads are on target" | sed -e 's/^[ \t]*//'

echo "-------------DEPTH/NUCLEOTIDES SEQUENCED INFO-------------"
python count_mutpos.py -mu $1.dcs.clipped.no_overlap.region.pileup -exmu $1.dcs.clipped.no_overlap.region.EXONS_ONLY.pileup -inmu $1.dcs.clipped.no_overlap.INTRONS_ONLY.region.pileup
