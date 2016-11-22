#!/bin/sh

echo "***PAIRED END UNIFIED CONSENSUS STATISTICS***"
echo " "

echo "-------------RAW READS-------------"
samtools flagstat $1.temp.sort.bam
echo " "

echo "-------------SSCS-------------"
samtools flagstat $1.sscs.sort.bam
echo " "

echo "-------------DCS-------------"
samtools flagstat $1.dcs.sort.bam
echo " "

echo "-------------DCS FILTERING-------------"
samtools flagstat $1.dcs.sort.mapped.q60.bam | awk 'NR==1{print $1 " Reads remaining after filtering for perfect mapping quality" ;}'
echo " "
samtools flagstat $1.dcs.filt.bam | awk 'NR==1{print $1 " Reads remaining after filtering for iSize" ;}'
