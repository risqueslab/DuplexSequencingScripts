#!/bin/bash

#	SSCS_PE_Unified_MITO.sh
#
#
#	Analogous to DS_PE_Unified_Mito Duplex Sequencing analysis pipeline using paired end reads and Scott Kennedy's Unified Consensus Maker,
# 	however only processing SSCS data for analysis. ***Uses an adjusted version of Count-Muts which account for read orientation.
#	
#	Modified by Dana Nachmanson 10/28/2016 to be used for the MITO analysis script.

clear

# Stop on any error inside or outside pipeline or on an unassigned variable.
set -e
set -o pipefail
set -u

# 1. SET RUN VARIABLES
tagLen=10          # Adaptor sequence length
spacerLen=1         # Spacer sequence length
readLen=100         # Sequencer read length
clipBegin=10         # Number of bases to clip off of beginning of reads
clipEnd=10           # Number of bases to clip off of end of reads

# 2. SET FILE LOCATIONS AND PATHS
DS_PATH=/Users/RRisques/Desktop/Duplex_Sequencing
PICARD_PATH=/Applications/Utilities/Seq_Analysis_Tools/picard-tools-2.2.1/picard.jar
ALIGN_REF=/Users/RRisques/Desktop/Duplex_Sequencing/Reference/Mito/MitoHumanAdjusted/human_g1k_v37.fasta
GATK_PATH=/Applications/Utilities/Seq_Analysis_Tools/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar
REGION_BED=/Users/RRisques/Desktop/Duplex_Sequencing/Reference/Mito/Human_mito.bed.txt # Bed file with genomic region
TRANS_TABLE=/Users/RRisques/Desktop/Duplex_Sequencing/Reference/Mito/VertebrateMito.transtable
MITO_REF=/Users/RRisques/Desktop/Duplex_Sequencing/Reference/Mito/MT2.fasta
MITO_REFERENCE=/Users/RRisques/Desktop/Duplex_Sequencing/Reference/Mito/MT2Reference.txt

# Adjusting variables for script use
finalReadLen=$((readLen-tagLen-spacerLen))
endTrimStart=$(($finalReadLen-$clipEnd+1))
iSize=$finalReadLen # Minimum insert size to allow through analysis
					# Purpose is to filter out small fragments that read into adaptors

# 3. SET SAMPLES TO ANALYZE:
# Folders need to contain the SSCS files created in original analysis, named SAMPLENAME.sscs.sort.bam and the SAMPLENAME.sscs.sort.bam.bai

folderList='C21 '

# 4. RUN SCRIPT:
# From the terminal-
# >> cd into the directory containing your SAMPLE FOLDERS and a copy of this script
# >> bash -x SSCS_PE_Unified_MITO.sh 2> SSCS_PE_Unified_Record_Mito.se  

for elmt in $folderList
do	
	cd ${elmt}
	
	# ***Currently using SSCS file aligned with BWA MEM***
	
	# Make a new folder named "SSCS_Analysis" in your sample folder and move aligned SSCS file there.
	mkdir SSCS_Analysis
	cp ${elmt}_mem.sscs.sort.bam ./SSCS_Analysis
	cp ${elmt}_mem.sscs.sort.bam.bai ./SSCS_Analysis
	
	cd SSCS_Analysis
	
    # Filter out unmapped reads
    samtools view -F4 ${elmt}_mem.sscs.sort.bam | samtools view -Sb -T $ALIGN_REF - > ${elmt}_mem.sscs.sort.mapped.bam
	samtools index ${elmt}_mem.sscs.sort.mapped.bam

	# Make file with only FORWARD mapping reads in a pair -------->
	samtools view -b -f 1 -F 16 ${elmt}_mem.sscs.sort.mapped.bam  > ${elmt}_mem.sscs.mapped.fwd.bam
	samtools index ${elmt}_mem.sscs.mapped.fwd.bam
	
	# Make file with only REVERSE mapping reads in a pair <--------
	samtools view -b -f 17 ${elmt}_mem.sscs.sort.mapped.bam  > ${elmt}_mem.sscs.mapped.rev.bam
	samtools index ${elmt}_mem.sscs.mapped.rev.bam

	# Make folder for each and move each plus index file in
	mkdir forward
	mkdir reverse
	mv ${elmt}_mem.sscs.mapped.fwd* ./forward 
	mv ${elmt}_mem.sscs.mapped.rev* ./reverse 

	# Process FORWARD mapping SSCS reads -------->
	# Go for forward directory
	cd forward
	
    # Filter for iSize
	samtools view ${elmt}_mem.sscs.mapped.fwd.bam | awk -v iSize="$iSize" -F'\t' '{if ($9 > iSize || $9 < -iSize) print $0}' > ${elmt}_mem.sscs.filt.sam
	samtools view -bT $ALIGN_REF ${elmt}_mem.sscs.filt.sam > ${elmt}_mem.sscs.filt.bam

    # End clipping SSCS reads
	java -jar -Xmx2g $PICARD_PATH AddOrReplaceReadGroups INPUT=${elmt}_mem.sscs.filt.bam OUTPUT=${elmt}.sscs.filt.readgroups.bam RGLB=Eeny RGPL=Meeny RGPU=Miny RGSM=Moe
	samtools index ${elmt}.sscs.filt.readgroups.bam
    java -jar -Xmx8g $GATK_PATH -T ClipReads -I ${elmt}.sscs.filt.readgroups.bam -o ${elmt}.sscs.filt.clipped.bam -R $ALIGN_REF --cyclesToTrim "$endTrimStart"-"$finalReadLen,1"-"$clipBegin" --clipRepresentation SOFTCLIP_BASES
	
    # Locally realign clipped SSCS reads, must be done for overlap clipping to work, --maxReadsForConsensus for realignment at high depth
    java -Xmx4g  -jar $GATK_PATH -T RealignerTargetCreator -dfrac 1 -R $ALIGN_REF -I ${elmt}.sscs.filt.clipped.bam -o ${elmt}.sscs.filt.clipped.bam.readgroups.intervals
    java -Xmx4g -jar $GATK_PATH -T IndelRealigner -dfrac 1 -R $ALIGN_REF -I ${elmt}.sscs.filt.clipped.bam -targetIntervals ${elmt}.sscs.filt.clipped.bam.readgroups.intervals -o ${elmt}.sscs.filt.clipped.realign.bam

    # Clip overlapping nucleotides in read pairs
    bam clipOverlap --in ${elmt}.sscs.filt.clipped.realign.bam --out ${elmt}.sscs.filt.clipped.realign.no_overlap.bam --stats 2> ${elmt}.overlapping-reads.stats

    # Create pileup. (Note: -Q0 filter is apparently critical for making mpileup work on softclipped reads)
    samtools mpileup -Q0 -B -A -d 500000 -f $ALIGN_REF ${elmt}.sscs.filt.clipped.realign.no_overlap.bam > ${elmt}.sscs.clipped.no_overlap.pileup

    # Generate final SSCS mutation data files:
    
    # Below depth -d is 50 and -C is 1 (turn off clonality cutoff). N filter 1 (no max). No -u (keep all, not just uniques). No minimum amount of mutations. All chromosomes.
    cat ${elmt}.sscs.clipped.no_overlap.pileup | python $DS_PATH/Programs/mut-position.1.31.py -C 1 -d 50 > ${elmt}.SSCS.all.for.mutpos

	# Filter out all positions that are not mito
	awk '{ if ($1 == "MT2" ) print }' ${elmt}.SSCS.all.for.mutpos > ${elmt}.SSCS.for.mutpos
    
    awk '{ if ($1 == "MT2" ) print }' ${elmt}.sscs.clipped.no_overlap.pileup > ${elmt}.sscs.clipped.no_overlap.onlyMITO.pileup
    
    # -C is 1 (turn off clonality cutoff). N filter 1 (no max). No -u (keep all, not just uniques). No minimum amount of mutations. JUST MITO
    cat ${elmt}.sscs.clipped.no_overlap.onlyMITO.pileup | python $DS_PATH/Programs/mut-position.1.31.py -C 1 > ${elmt}.SSCS.for.mutpos
    
	# Generate count_muts file. -C is 1 (turn off clonality cutoff). No -u (keep all, not just uniques). No minimum amount of mutations.
    cat ${elmt}.sscs.clipped.no_overlap.onlyMITO.pileup | python $DS_PATH/Programs/CountMuts.1.41.py -C 1 > ${elmt}.SSCS.for.countmuts
    
    # Mito mutpos post processing
	python $DS_PATH/Programs/Mito/MutationConsequences2.py $TRANS_TABLE $REGION_BED $MITO_REF ${elmt}.SSCS.for.mutpos ${elmt} 0 1
	python $DS_PATH/Programs/Mito/mutpos_analyzer.py ${elmt}.SSCS.for.mutpos ${elmt}.consequences.txt ${elmt} $MITO_REFERENCE
	
	cd ..
	
	# Process REVERSE mapping SSCS reads -------->
	# Go for reverse directory
	cd reverse
	
    # Filter for iSize
	samtools view ${elmt}_mem.sscs.mapped.rev.bam | awk -v iSize="$iSize" -F'\t' '{if ($9 > iSize || $9 < -iSize) print $0}' > ${elmt}_mem.sscs.filt.sam
	samtools view -bT $ALIGN_REF ${elmt}_mem.sscs.filt.sam > ${elmt}_mem.sscs.filt.bam

    # End clipping SSCS reads
	java -jar -Xmx2g $PICARD_PATH AddOrReplaceReadGroups INPUT=${elmt}_mem.sscs.filt.bam OUTPUT=${elmt}.sscs.filt.readgroups.bam RGLB=Eeny RGPL=Meeny RGPU=Miny RGSM=Moe
	samtools index ${elmt}.sscs.filt.readgroups.bam
    java -jar -Xmx8g $GATK_PATH -T ClipReads -I ${elmt}.sscs.filt.readgroups.bam -o ${elmt}.sscs.filt.clipped.bam -R $ALIGN_REF --cyclesToTrim "$endTrimStart"-"$finalReadLen,1"-"$clipBegin" --clipRepresentation SOFTCLIP_BASES
	
    # Locally realign clipped SSCS reads, must be done for overlap clipping to work, --maxReadsForConsensus for realignment at high depth
    java -Xmx4g  -jar $GATK_PATH -T RealignerTargetCreator -dfrac 1 -R $ALIGN_REF -I ${elmt}.sscs.filt.clipped.bam -o ${elmt}.sscs.filt.clipped.bam.readgroups.intervals
    java -Xmx4g -jar $GATK_PATH -T IndelRealigner -dfrac 1 -R $ALIGN_REF -I ${elmt}.sscs.filt.clipped.bam -targetIntervals ${elmt}.sscs.filt.clipped.bam.readgroups.intervals -o ${elmt}.sscs.filt.clipped.realign.bam

    # Clip overlapping nucleotides in read pairs
    bam clipOverlap --in ${elmt}.sscs.filt.clipped.realign.bam --out ${elmt}.sscs.filt.clipped.realign.no_overlap.bam --stats 2> ${elmt}.overlapping-reads.stats

    # Create pileup. (Note: -Q0 filter is apparently critical for making mpileup work on softclipped reads)
    samtools mpileup -Q0 -B -A -d 500000 -f $ALIGN_REF ${elmt}.sscs.filt.clipped.realign.no_overlap.bam > ${elmt}.sscs.clipped.no_overlap.pileup

    # Generate final SSCS mutation data files:
    
    # -C is 1 (turn off clonality cutoff). N filter 1 (no max). No -u (keep all, not just uniques). No minimum amount of mutations. No minimum depth. All chromosomes.
    # Use an adjusted version of mut-position to adjust for reverse complemented bases.
    cat ${elmt}.sscs.clipped.no_overlap.pileup | python $DS_PATH/Programs/SSCS_Analysis/mut-position.1.31_antiref.py -d 1 -C 1 > ${elmt}.SSCS.all.rev.mutpos

	# Filter out all positions that are not mito
	awk '{ if ($1 == "MT2" ) print }' ${elmt}.SSCS.all.rev.mutpos > ${elmt}.SSCS.rev.mutpos
    
    awk '{ if ($1 == "MT2" ) print }' ${elmt}.sscs.clipped.no_overlap.pileup > ${elmt}.sscs.clipped.no_overlap.onlyMITO.pileup
    
    # -C is 1 (turn off clonality cutoff). N filter 1 (no max). No -u (keep all, not just uniques). No minimum amount of mutations. JUST MITO. No minimum depth.
    cat ${elmt}.sscs.clipped.no_overlap.onlyMITO.pileup | $DS_PATH/Programs/SSCS_Analysis/mut-position.1.31_antiref.py -d 1 -C 1 > ${elmt}.SSCS.rev.mutpos
    
	# Generate count_muts file. -C is 1 (turn off clonality cutoff). No -u (keep all, not just uniques). No minimum amount of mutations. No minimum depth.
	# Use an adjusted version of countMuts to adjust for reverse complemented bases.
    cat ${elmt}.sscs.clipped.no_overlap.onlyMITO.pileup | python $DS_PATH/Programs/SSCS_Analysis/CountMuts.1.41_antiref.py -d 1 -C 1 > ${elmt}.SSCS.rev.countmuts
    
    # Mito mutpos post processing
	python $DS_PATH/Programs/Mito/MutationConsequences2.py $TRANS_TABLE $REGION_BED $MITO_REF ${elmt}.SSCS.rev.mutpos ${elmt} 0 1
	python $DS_PATH/Programs/Mito/mutpos_analyzer.py ${elmt}.SSCS.rev.mutpos ${elmt}.consequences.txt ${elmt} $MITO_REFERENCE
	
	cd ..
	cd ..
	cd ..
done