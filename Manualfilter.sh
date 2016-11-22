#!/bin/bash

# 1. SET RUN VARIABLES
tagLen=10           # Adapter sequence length
spacerLen=1         # Spacer sequence length
readLen=100         # Sequencer read length
clipBegin=10        # Number of bases to clip off of beginning of reads
clipEnd=10          # Number of bases to clip off of end of reads

# Adjusting variables for script use
finalReadLen=$((readLen-tagLen-spacerLen))
endTrimStart=$(($finalReadLen-$clipEnd+1))
iSize=$finalReadLen # Minimum insert size to allow through analysis
					# Purpose is to filter out small fragments that read into adaptors

# Filter for mapping only
samtools view -F4 B25.dcs.sort.bam | samtools view -Sb -T /Users/RRisques/Desktop/Duplex_Sequencing/Reference/Mito/MitoHumanAdjusted/human_g1k_v37.fasta - > MAPPING.bam
samtools index MAPPING.bam


# Filter out reads with less than perfect mapping scores (<60)
samtools view -bh -q60 MAPPING.bam > MAPPING.Q60.bam
samtools index MAPPING.Q60.bam


# Filter for size
samtools view MAPPING.Q60.bam | awk -v iSize="$iSize" -F'\t' '{if ($9 > iSize || $9 < -iSize) print $0}' > MAPPING.Q60.ISIZE.sam
