#!/bin/bash

#	SSCS_PE_Unified.sh
#
#
#	Analogous to DS_PE_Unified Duplex Sequencing analysis pipeline using paired end reads and Scott Kennedy's Unified Consensus Maker,
# 	however only processing SSCS data for analysis. ***Uses an adjusted version of Count-Muts which account for read orientation.
#	
#	Last modified by Dana Nachmanson 9/12/16.  

clear

# Stop on any error inside or outside pipeline or on an unassigned variable.
set -e
set -o pipefail
set -u

# 1. SET RUN VARIABLES
tagLen=12           # Adaptor sequence length
spacerLen=5         # Spacer sequence length
readLen=100         # Sequencer read length
clipBegin=5         # Number of bases to clip off of beginning of reads
clipEnd=5           # Number of bases to clip off of end of reads

# 2. SET FILE LOCATIONS AND PATHS
DS_PATH=/Users/RRisques/Desktop/Duplex_Sequencing/Programs
PICARD_PATH=/Applications/Utilities/Seq_Analysis_Tools/picard-tools-2.2.1/picard.jar
ALIGN_REF=/Users/RRisques/Desktop/Duplex_Sequencing/Reference/Human_Genome/human_g1k_v37.fasta
GATK_PATH=/Applications/Utilities/Seq_Analysis_Tools/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar
REGION_BED=/Users/RRisques/Desktop/Duplex_Sequencing/Reference/P53/p53region.bed # Bed file with genomic region
SNP_BED=/Users/RRisques/Desktop/Duplex_Sequencing/Reference/P53/p53SNPs.bed      # Bed file with SNP positions 

# Adjusting variables for script use
finalReadLen=$((readLen-tagLen-spacerLen))
endTrimStart=$(($finalReadLen-$clipEnd+1))
iSize=$finalReadLen # Minimum insert size to allow through analysis
					# Purpose is to filter out small fragments that read into adaptors

# 3. SET SAMPLES TO ANALYZE:
# Folders need to contain the SSCS files created in original analysis, named SAMPLENAME.sscs.sort.bam and the SAMPLENAME.sscs.sort.bam.bai

folderList='GRE'

# 4. RUN SCRIPT:
# From the terminal-
# >> cd into the directory containing your SAMPLE FOLDERS and a copy of this script
# >> bash -x SSCS_PE_Unified.sh 2> SSCS_PE_Unified_Record.se  

for elmt in $folderList
do	
	cd ${elmt}
	
	# ***Currently using SSCS file aligned with BWA SW and not touching MEM file***
	
	# Make a new folder named "SSCS_Analysis" in your sample folder and move aligned SSCS file there.
	mkdir SSCS_Analysis
	cp *.sscs.sort.bam ./SSCS_Analysis
	cp *.sscs.sort.bai ./SSCS_Analysis
	
	cd SSCS_Analysis

    # Filter out unmapped reads
    samtools view -F4 ${elmt}.sscs.sort.bam | samtools view -Sb -T $ALIGN_REF - > ${elmt}.sscs.sort.mapped.bam
	samtools index ${elmt}.sscs.sort.mapped.bam
	
	# Filter out reads with less than perfect mapping scores (<60)
	samtools view -bh -q60 ${elmt}.sscs.sort.mapped.bam > ${elmt}.sscs.sort.mapped.q60.bam
	samtools index ${elmt}.sscs.sort.mapped.q60.bam
    
    # Filter for size
	samtools view ${elmt}.sscs.sort.mapped.q60.bam | awk -v iSize="$iSize" -F'\t' '{if ($9 > iSize || $9 < -iSize) print $0}' > ${elmt}.sscs.filt.sam
	samtools view -bT $ALIGN_REF ${elmt}.sscs.filt.sam > ${elmt}.sscs.filt.bam

    # End clipping SSCS reads
	java -jar -Xmx2g $PICARD_PATH AddOrReplaceReadGroups INPUT=${elmt}.sscs.filt.bam OUTPUT=${elmt}.sscs.filt.readgroups.bam RGLB=Eeny RGPL=Meeny RGPU=Miny RGSM=Moe
	samtools index ${elmt}.sscs.filt.readgroups.bam
    java -jar -Xmx8g $GATK_PATH -T ClipReads -I ${elmt}.sscs.filt.readgroups.bam -o ${elmt}.sscs.filt.clipped.bam -R $ALIGN_REF --cyclesToTrim "$endTrimStart"-"$finalReadLen,1"-"$clipBegin" --clipRepresentation SOFTCLIP_BASES
	
    # Locally realign clipped SSCS reads, must be done for overlap clipping to work, --maxReadsForConsensus for realignment at high depth
    java -Xmx4g  -jar $GATK_PATH -T RealignerTargetCreator -dfrac 1 -R $ALIGN_REF -I ${elmt}.sscs.filt.clipped.bam -o ${elmt}.sscs.filt.clipped.bam.readgroups.intervals
    java -Xmx4g -jar $GATK_PATH -T IndelRealigner -dfrac 1 -R $ALIGN_REF -I ${elmt}.sscs.filt.clipped.bam -targetIntervals ${elmt}.sscs.filt.clipped.bam.readgroups.intervals -o ${elmt}.sscs.filt.clipped.realign.bam

    # Clip overlapping nucleotides in read pairs
    bam clipOverlap --in ${elmt}.sscs.filt.clipped.realign.bam --out ${elmt}.sscs.filt.clipped.realign.no_overlap.bam --stats 2> ${elmt}.overlapping-reads.stats

    # Create pileup. (Note: -Q0 filter is apparently critical for making mpileup work on softclipped reads)
    samtools mpileup -Q0 -B -A -d 500000 -f $ALIGN_REF ${elmt}.sscs.filt.clipped.realign.no_overlap.bam > ${elmt}.sscs.clipped.no_overlap.pileup
    
    # Create pileup with JUST genomic coordinates of interest
    python $DS_PATH/filter_pileup.py $REGION_BED ${elmt}.sscs.clipped.no_overlap.pileup ${elmt}.sscs.clipped.no_overlap.region.pileup

    # Create pileup file without SNPs
    awk 'FNR==NR { a[$1 '\t' $2]; next } !("chr" $1 '\t' $2 in a)' $SNP_BED ${elmt}.sscs.clipped.no_overlap.region.pileup > ${elmt}.sscs.clipped.no_overlap.region.noSNPs.pileup
    
    # Generate final SSCS mutation data file:
    
	# Generate count_muts file. Below depth -d is 1 (count all sites regardless of depth) and -C is 1 (turn off clonality cutoff). Default N filter is 5%. -u is counting uniques only .
    cat ${elmt}.sscs.clipped.no_overlap.region.noSNPs.pileup | python $DS_PATH/countMuts.1.41_SSCS.py -d 1 -C 1 -u > ${elmt}.SSCS.noSNPs.pileup.countmuts
    
    # Generate muts file, minimum depth is 1, minimum number of mutations is 1. No clonality filter by default.
    cat ${elmt}.sscs.clipped.no_overlap.region.noSNPs.pileup | python $DS_PATH/mut-position.1.31.py -d 1 -n 1 > ${elmt}.sscs-muts.noSNP.txt
	
	# Generate mutpos file, minimum depth is 1, minimum number of mutations is 0. No clonality filter by default.
	cat ${elmt}.sscs.clipped.no_overlap.region.noSNPs.pileup | python $DS_PATH/mut-position.1.31.py -d 1 > ${elmt}.sscs_clipped.bam.pileup.noSNP.mutpos

	
	cd ..
	cd ..
done