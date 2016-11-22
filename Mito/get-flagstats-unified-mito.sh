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
samtools flagstat $1.dcs.filt.bam | awk 'NR==1{print $1 " Reads remaining after filtering for iSize" ;}'
echo "$(samtools view $1.dcs.filt.clipped.realign.bam "MT2" | wc -l) Reads are on target" | sed -e 's/^[ \t]*//'
