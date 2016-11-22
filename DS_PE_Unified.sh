#!/bin/bash

#	DS_PE_Unified.sh
#
#
#	Duplex Sequencing analysis pipeline using paired end reads and Scott Kennedy's Unified Consensus Maker.
# 	
#	
#	Last modified by Dana Nachmanson 9/2/16: 
#		-Removed muts by position plot and replaced with mutations per cycle plot made with GATK
#		-Added a depth by position plot
#		-Added histogram of iSizes
#		-Set iSize to be a value decided by read length and tag+spacer length 
#		-Added an mutpos output file that spits out only SNPs with no min or max # of muts
#		-Replaced downsampling flagging on GATK commands
#		-Removed -fix_misencoded_quality scores flag from GATK 
#
#

clear

# Stop on any error inside or outside pipeline or on an unassigned variable.
set -e
set -o pipefail
set -u

# 1. SET RUN VARIABLES
minMem=3            # Minimum number of reads to reach consensus
maxMem=200         # Maximum number of reads to reach consensus
cutOff=0.7          # % of nucleotides at a position in read that must be identical in order for consensus at that position
nCutOff=1           # Maximum fraction of Ns allowed in a consensus
tagLen=10           # Adapter sequence length
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
# Folder names separated by spaces containing R1 and R2 Fastq.gz files named: samplename.seq1.fastq.gz and samplename.seq2.fastq.gz
# NOTE: This script can use compressed Fastq files, no need to unzip files before running this script

folderList='GRE'

# 4. RUN SCRIPT:
# From the terminal-
# >> cd into the directory containing your SAMPLE FOLDERS and a copy of this script
# >> bash -x DS_PE_Unified.sh 2> DS_PE_Unified_Record.se  

for elmt in $folderList
do	
	cd ${elmt}
	# Consensus Maker
    java -jar $PICARD_PATH FastqToSam F1=${elmt}.seq1.fastq.gz F2=${elmt}.seq2.fastq.gz O=/dev/stdout SM=${elmt}|python $DS_PATH/UnifiedConsensusMaker.py --input /dev/stdin --taglen ${tagLen} --spacerlen ${spacerLen} --minmem ${minMem} --maxmem ${maxMem} --cutoff ${cutOff} --Ncutoff ${nCutOff} --write-sscs --prefix ${elmt} --tagstats

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

    #***Currently using DCS files aligned with BWA MEM only ***

    # Filter out unmapped reads
    samtools view -F4 ${elmt}_mem.dcs.sort.bam | samtools view -Sb -T $ALIGN_REF - > ${elmt}.dcs.sort.mapped.bam
	samtools index ${elmt}.dcs.sort.mapped.bam
	
	# Filter out reads with less than perfect mapping scores (<60)
	samtools view -bh -q60 ${elmt}.dcs.sort.mapped.bam > ${elmt}.dcs.sort.mapped.q60.bam
	samtools index ${elmt}.dcs.sort.mapped.q60.bam
    
    # Filter for insert size
	samtools view ${elmt}.dcs.sort.mapped.q60.bam | awk -v iSize="$iSize" -F'\t' '{if ($9 > iSize || $9 < -iSize) print $0}' > ${elmt}.dcs.filt.sam
	samtools view -bT $ALIGN_REF ${elmt}.dcs.filt.sam > ${elmt}.dcs.filt.bam

    # End clipping DCS reads
	java -jar -Xmx2g $PICARD_PATH AddOrReplaceReadGroups INPUT=${elmt}.dcs.filt.bam OUTPUT=${elmt}.dcs.filt.readgroups.bam RGLB=Eeny RGPL=Meeny RGPU=Miny RGSM=Moe
	samtools index ${elmt}.dcs.filt.readgroups.bam
    java -jar -Xmx8g $GATK_PATH -T ClipReads -I ${elmt}.dcs.filt.readgroups.bam -o ${elmt}.dcs.filt.clipped.bam -R $ALIGN_REF --cyclesToTrim "$endTrimStart"-"$finalReadLen,1"-"$clipBegin" --clipRepresentation SOFTCLIP_BASES
	
    # Locally realign clipped DCS reads, must be done for overlap clipping to work, --maxReadsForConsensus for realignment at high depth
    java -Xmx4g  -jar $GATK_PATH -T RealignerTargetCreator -dfrac 1 -R $ALIGN_REF -I ${elmt}.dcs.filt.clipped.bam -o ${elmt}.dcs.filt.clipped.bam.readgroups.intervals
    java -Xmx4g -jar $GATK_PATH -T IndelRealigner -dfrac 1 -R $ALIGN_REF -I ${elmt}.dcs.filt.clipped.bam -targetIntervals ${elmt}.dcs.filt.clipped.bam.readgroups.intervals -o ${elmt}.dcs.filt.clipped.realign.bam

    # Clip overlapping nucleotides in read pairs
    bam clipOverlap --in ${elmt}.dcs.filt.clipped.realign.bam --out ${elmt}.dcs.filt.clipped.realign.no_overlap.bam --stats 2> ${elmt}.overlapping-reads.stats

    # Create pileup. (Note: -Q0 filter is apparently critical for making mpileup work on softclipped reads)
    samtools mpileup -Q0 -B -A -d 500000 -f $ALIGN_REF ${elmt}.dcs.filt.clipped.realign.no_overlap.bam > ${elmt}.dcs.clipped.no_overlap.pileup
    
    # Create pileup with JUST genomic coordinates of interest
    python $DS_PATH/filter_pileup.py $REGION_BED ${elmt}.dcs.clipped.no_overlap.pileup ${elmt}.dcs.clipped.no_overlap.region.pileup N

    # Create pileup file without SNPs
    awk 'FNR==NR { a[$1 '\t' $2]; next } !("chr" $1 '\t' $2 in a)' $SNP_BED ${elmt}.dcs.clipped.no_overlap.region.pileup > ${elmt}.dcs.clipped.no_overlap.region.noSNPs.pileup
    
    # Create pileup file with ONLY SNPs
    awk 'FNR==NR { a[$1 '\t' $2]; next } ("chr" $1 '\t' $2 in a)' $SNP_BED ${elmt}.dcs.clipped.no_overlap.region.pileup > ${elmt}.dcs.clipped.no_overlap.region.ONLYSNPs.pileup

    
   # Generate final DCS mutation data files:
    
    # Final files with SNPs:
    # Minimum depth of 1, minimum number of mutations is 1. By default, no clonality filters.
    cat ${elmt}.dcs.clipped.no_overlap.region.pileup | python $DS_PATH/mut-position.1.31.py -d 1 -n 1 > ${elmt}.DCS-muts.txt
    
    # Minimum depth of 1. By default, minimum number of mutations is 0 and no clonality filters.
    cat ${elmt}.dcs.clipped.no_overlap.region.pileup | python $DS_PATH/mut-position.1.31.py -d 1 > ${elmt}.DCS.pileup.mutpos
	
	# Final files without SNPs:
	# Below depth -d is 1 (count all sites regardless of depth) and -C is 1 (turn off clonality cutoff). Default N filter is 5%. -u is counting uniques only .
    cat ${elmt}.dcs.clipped.no_overlap.region.noSNPs.pileup | python $DS_PATH/countMuts.py -d 1 -C 1 -u > ${elmt}.DCS.noSNPs.pileup.countmuts
	
	# Minimum depth of 1, minimum number of mutations is 1. By default, no clonality filters.
    cat ${elmt}.dcs.clipped.no_overlap.region.noSNPs.pileup | python $DS_PATH/mut-position.1.31.py -d 1 -n 1 > ${elmt}.DCS-muts.noSNPs.txt
    
    # Minimum depth of 1. By default, minimum number of mutations is 0 and no clonality filters.
    cat ${elmt}.dcs.clipped.no_overlap.region.noSNPs.pileup | python $DS_PATH/mut-position.1.31.py -d 1 > ${elmt}.DCS.noSNPs.pileup.mutpos
	
	# Final files with ONLY SNPs:
	# Minimum depth of 1. By default, minimum number of mutations is 0 and no clonality filters.
    cat ${elmt}.dcs.clipped.no_overlap.region.ONLYSNPs.pileup | python $DS_PATH/mut-position.1.31.py -d 1 > ${elmt}.DCS.ONLYSNPs.pileup.mutpos

    # Generate statistics files:
    
    # Plot DCS depth by genomic coordinate
    python $DS_PATH/Plot/plot_depth_by_position.py ${elmt}.dcs.clipped.no_overlap.region.pileup
    
    # Plot insert-size histogram (using unfiltered and unclipped data)
	java -jar $PICARD_PATH CollectInsertSizeMetrics I=${elmt}_mem.dcs.sort.bam O=${elmt}.iSize_Metrics.txt H=${elmt}.iSize_Histogram.pdf M=0.5
	
	# Plot mutations by read cycle 
	java -jar /$GATK_PATH -T ErrorRatePerCycle -R $ALIGN_REF -dfrac 1 -I ${elmt}.dcs.filt.readgroups.bam -o ${elmt}.ErrorRatePerCycle.txt 
	python $DS_PATH/Plot/plot_error_by_cycle.py ${elmt}.ErrorRatePerCycle.txt
	
    # Print reads statistics
    bash $DS_PATH/get-flagstats-unified.sh ${elmt} > ${elmt}.flagstats.stats

	cd ..
done