#!/bin/bash

#	DS_PE_Unified_Mito.sh
#
#
#	Duplex Sequencing analysis pipeline using paired end reads and a unified duplex consensus maker.
#	8/12/16 - Dana Nachmanson adjusted script to customize to mito DNA analysis pipeline 
# 	10/24/16 - Dana Nachmanson removed JSON filter and added explicit filtering of BAM files

# Stop on any error inside or outside pipeline or on an unassigned variable.
set -e
set -o pipefail
set -u

# 1. SET RUN VARIABLES
minMem=3            # Minimum number of reads to reach consensus
maxMem=1000         # Maximum number of reads to reach consesnsus
cutOff=0.7          # % of nucleotides at a position in read that must be identical in order for consensus at that position
nCutOff=0.05        # Maximum fraction of Ns allowed in a consensus
tagLen=10           # Adaptor sequence length
spacerLen=1         # Spacer sequence length
readLen=100         # Sequencer read length
clipBegin=10        # Number of bases to clip off of beginning of reads
clipEnd=10          # Number of bases to clip off of end of reads

# 2. SET FILE LOCATIONS AND PATHS
DS_PATH=/Users/RRisques/Desktop/Duplex_Sequencing
PICARD_PATH=/Applications/Utilities/Seq_Analysis_Tools/picard-tools-2.2.1/picard.jar
ALIGN_REF=/Users/RRisques/Desktop/Duplex_Sequencing/Reference/Mito/MitoHumanAdjusted/human_g1k_v37.fasta
GATK_PATH=/Applications/Utilities/Seq_Analysis_Tools/GenomeAnalysisTK-3.6
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
# Folder names separated by spaces containing R1 and R2 Fastq.gz files named: samplename.seq1.fastq.gz and samplename.seq2.fastq.gz
# NOTE: This script can use compressed Fastq files, no need to unzip files before running this script

folderList='5916GAY'

# 4. RUN SCRIPT:
# From the terminal-
# >> cd into the directory containing your SAMPLE FOLDERS and a copy of this script
# >> bash -x DS_PE_Unified_Mito.sh 2> DS_PE_Unified_Mito_Record.se  

for elmt in $folderList
do	
	cd ${elmt}
	# Consensus Maker
    java -jar $PICARD_PATH FastqToSam F1=${elmt}.seq1.fastq.gz F2=${elmt}.seq2.fastq.gz O=/dev/stdout SM=${elmt}|python $DS_PATH/Programs/UnifiedConsensusMaker.py --input /dev/stdin --taglen ${tagLen} --spacerlen ${spacerLen} --minmem ${minMem} --maxmem ${maxMem} --cutoff ${cutOff} --Ncutoff ${nCutOff} --write-sscs --prefix ${elmt} --tagstats

    # Align forward and reverse SSCS reads using both BWA algorithms SW and MEM
    bwa aln $ALIGN_REF ${elmt}_read1_sscs.fq.gz > ${elmt}.sscs.read1.aln
    bwa aln $ALIGN_REF ${elmt}_read2_sscs.fq.gz > ${elmt}.sscs.read2.aln

    bwa sampe -s $ALIGN_REF ${elmt}.sscs.read1.aln ${elmt}.sscs.read2.aln ${elmt}_read1_sscs.fq.gz ${elmt}_read2_sscs.fq.gz | samtools sort -o ${elmt}.sscs.sort.bam -

    bwa mem $ALIGN_REF ${elmt}_read1_sscs.fq.gz ${elmt}_read2_sscs.fq.gz|samtools sort -o ${elmt}_mem.sscs.sort.bam -

	# Align forward and reverse DCS reads using both BWA algorithms SW and MEM
    bwa aln $ALIGN_REF ${elmt}_read1_dcs.fq.gz > ${elmt}.dcs.read1.aln
    bwa aln $ALIGN_REF ${elmt}_read2_dcs.fq.gz > ${elmt}.dcs.read2.aln
    bwa sampe -s $ALIGN_REF ${elmt}.dcs.read1.aln ${elmt}.dcs.read2.aln ${elmt}_read1_dcs.fq.gz ${elmt}_read2_dcs.fq.gz | samtools sort -o ${elmt}.dcs.sort.bam -
    bwa mem $ALIGN_REF ${elmt}_read1_dcs.fq.gz ${elmt}_read2_dcs.fq.gz|samtools sort -o ${elmt}_mem.dcs.sort.bam -

    # Index both SSCS and DCS files
	samtools index ${elmt}.sscs.sort.bam
    samtools index ${elmt}_mem.sscs.sort.bam
    samtools index ${elmt}.dcs.sort.bam
    samtools index ${elmt}_mem.dcs.sort.bam

    #***Currently using DCS aligned with MEM file and NO mapping score filter***

 	# Filter out unmapped reads
    samtools view -F4 ${elmt}_mem.dcs.sort.bam | samtools view -Sb -T $ALIGN_REF - > ${elmt}_mem.dcs.sort.mapped.bam
	samtools index ${elmt}_mem.dcs.sort.mapped.bam
    
    # Filter for insert size
	samtools view ${elmt}_mem.dcs.sort.mapped.bam | awk -v iSize="$iSize" -F'\t' '{if ($9 > iSize || $9 < -iSize) print $0}' > ${elmt}_mem.dcs.filt.sam
	samtools view -bT $ALIGN_REF ${elmt}_mem.dcs.filt.sam > ${elmt}.dcs.filt.bam

    # End clipping DCS reads
	java -jar -Xmx2g $PICARD_PATH AddOrReplaceReadGroups INPUT=${elmt}.dcs.filt.bam OUTPUT=${elmt}.dcs.filt.readgroups.bam RGLB=Eeny RGPL=Meeny RGPU=Miny RGSM=Moe
	samtools index ${elmt}.dcs.filt.readgroups.bam
    java -jar -Xmx8g $GATK_PATH/GenomeAnalysisTK.jar -T ClipReads -I ${elmt}.dcs.filt.readgroups.bam -o ${elmt}.dcs.filt.clipped.bam -R $ALIGN_REF --cyclesToTrim "$endTrimStart"-"$finalReadLen,1"-"$clipBegin" --clipRepresentation SOFTCLIP_BASES --fix_misencoded_quality_scores
	
    # Locally realign clipped DCS reads, must be done for overlap clipping to work, --maxReadsForConsensus for realignment at high depth
    java -Xmx4g  -jar $GATK_PATH/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $ALIGN_REF -I ${elmt}.dcs.filt.clipped.bam -o ${elmt}.dcs.filt.clipped.bam.readgroups.intervals
    java -Xmx4g -jar $GATK_PATH/GenomeAnalysisTK.jar -T IndelRealigner --maxReadsForConsensuses 100000000 -R $ALIGN_REF -I ${elmt}.dcs.filt.clipped.bam -targetIntervals ${elmt}.dcs.filt.clipped.bam.readgroups.intervals -o ${elmt}.dcs.filt.clipped.realign.bam

    # Clip overlapping nucleotides in read pairs
    bam clipOverlap --in ${elmt}.dcs.filt.clipped.realign.bam --out ${elmt}.dcs.filt.clipped.realign.no_overlap.bam --stats 2> ${elmt}.overlapping-reads.stats
    
    # Create pileup. (Note: -Q0 filter is apparently critical for making mpileup work on softclipped reads)
    samtools mpileup -Q0 -B -A -d 500000 -f $ALIGN_REF ${elmt}.dcs.filt.clipped.realign.no_overlap.bam > ${elmt}.dcs.clipped.no_overlap.pileup 
    
    # Generate final DCS mutation data files:
    
    # Below depth -d is 50 and -C is 1 (turn off clonality cutoff). N filter 1 (no max). No -u (keep all, not just uniques). No minimim amount of mutations.
    cat ${elmt}.dcs.clipped.no_overlap.pileup | python $DS_PATH/Programs/mut-position.1.31.py -C 1 -d 50 > ${elmt}.DCS.all.mutpos

	# Filter out all positions that are not mito
	awk '{ if ($1 == "MT2" ) print }' ${elmt}.DCS.all.mutpos  > ${elmt}.DCS.mutpos
	
	# Mito mutpos post processing
	python $DS_PATH/Programs/Mito/MutationConsequences2.py $TRANS_TABLE $REGION_BED $MITO_REF ${elmt}.DCS.mutpos ${elmt} 0 1
	python $DS_PATH/Programs/Mito/mutpos_analyzer.py ${elmt}.DCS.mutpos ${elmt}.consequences.txt ${elmt} $MITO_REFERENCE
	
    # Generate statistics files:
    
    # Print reads
    bash $DS_PATH/Programs/Mito/get-flagstats-unified-mito.sh ${elmt} > ${elmt}.flagstats.txt ${elmt}
    
    # Plot depth by position using pileup filtered for only chromosome MT2
    awk '{ if ($1 == "MT2" ) print }' ${elmt}.dcs.clipped.no_overlap.pileup  > ${elmt}.dcs.clipped.only_MT2.pileup
	python $DS_PATH/Programs/Plot/plot_depth_by_position.py ${elmt}.dcs.clipped.only_MT2.pileup
    
    # Plot insert-size histogram (using unfiltered and unclipped data)
	java -jar $PICARD_PATH CollectInsertSizeMetrics I=${elmt}_mem.dcs.sort.bam O=${elmt}.iSize_Metrics.txt H=${elmt}.iSize_Histogram.pdf M=0.5

	cd ..
done