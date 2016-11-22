#!/bin/bash

# Generated on Tue Dec 31 2015 with the command by Kate Bayliss:

# 3/17-20/16 Many adjustments by J. Salk, statistics files, SNP filtering functionality
# 3/24/16 PM added back local realignment to fix overlap clipping
# 7/13/16 Added in SNP reporting function

# TO DO:
# Adapt for unified consensus maker
# For local re-alignment step need to add option to allow higher depth (see http://gatkforums.broadinstitute.org/gatk/discussion/3292/indelrealigner-not-attempting-realignment-in-interval)
# Need to add an iSize filtration step after alignment step

#-------------------------------------------------------------------------------------------------

# Step 1: Setup variables for run:
clear

# Stop on any error inside or outside pipeline or on an unassigned variable.
set -e
set -o pipefail
set -u

#DEFAULTS
DSpath=''
alignRef=''
runIdentifier=''
read1in=seq1.fq
read2in=seq2.fq
iSize=-1		### Filtering for iSize. -1 turns off
minMem=3
maxMem=1000		### once reaches 1000 keeps going but only retains a random selection of 1000 per family
cutOff=0.7
nCutOff=1		
readLength=0		
barcodeLength=12
spacerLength=5
filtersSet='os' 	### o, s, n are options. o = overlap filter (dumps read pair), s = softclipped out, n = if exceed N's
readTypes='dpm'		### see readme documentation. Selects for only certain types of flags (2nd colum in BAM file) 
repFilt=9
readOut=1000000

#NONDEFAULTS
DSpath='/Users/jjsalk/Dropbox/DS-software/Duplex-PE-2.1/programs'
PICARDPATH='/Users/jjsalk/Dropbox/DS-software/programs/picard-tools-1'
GATKPATH='/Users/jjsalk/Dropbox/DS-software/programs/GenomeAnalysisTK-3.4-46' 
PYTHONPATH=/Users/jjsalk/Dropbox/DS-software/programs/biopython-1.62b/			### Critical for Biopython modules to work

runIdentifier=seq
readLength=101		### 101 is for King lab sequencer
barcodeLength=12 	
spacerLength=5		
filtersSet=s   		### only filtering out softclipped reads; overlapping (o) and hamming filters are retained
readTypes=dm   		### only making consensus sequences from aligned reads. OK if mate is not aligned.

# --------------------------------------------------------------------------------------------------------------------------------

#SET THESE PARAMETERS:

alignRef='/Users/jjsalk/Dropbox/DS-software/refgenomes/human_g1k_v37/human_g1k_v37.fa'	
trimBeginning=12
trimEnd=12

# This trims pileup file to only exact regions of interest (i.e. exons, etc). To adjust need modify the called script)
TARGETREGION='/Users/jjsalk/Dropbox/DS-software/pileup-processing-scripts/filter_regions/pileup-target-filter-p53AR.script'

# First is SNP sites to exclude, 2nd set is SNP sites to keep (typically will be the same but coded differently)
SNPPURGE="/Users/jjsalk/Dropbox/DS-software/pileup-processing-scripts/filter_SNPs/p53AR/purge-p53AR-MAF-0.001.script"; export SNPPURGE
SNPSIFT="/Users/jjsalk/Dropbox/DS-software/pileup-processing-scripts/filter_SNPs/p53AR/sift-p53AR-MAF-0.001.script"; export SNPSIFT


# --------------------------------------------------------------------------------------------------------------------------------

#FINAL READ LENGTH
readLength=$((readLength-barcodeLength-spacerLength))

#TRIM POSITIONS
trimBeginningB=1
trimBeginningE=$trimBeginning
trimEndB=$(($readLength-$trimEnd+1))
trimEndE=$readLength

#LOG_FILE_NAME
logFile=${runIdentifier}.log.txt

#Export all variables
export DSpath
export PICARDPATH
export GATKPATH
export PYTHONPATH
export alignRef
export runIdentifier
export read1in
export read2in
export iSize
export minMem
export maxMem
export cutOff
export nCutOff
export readLength
export barcodeLength
export spacerLength
export filtersSet
export readTypes
export repFilt
export readOut
export trimBeginningB
export trimBeginningE
export trimEndB
export trimEndE
export TARGETREGION
export SNPPURGE
export SNPSIFT


# If put the two lines below around a block of code, it is equivilant to commenting out every line. Use to run select parts of script.
: <<'COMMENT'
COMMENT

# Print out options used to log file
touch $logFile
echo "Run identifier: " $runIdentifier | tee -a ${logFile}
echo "Program path: " $DSpath | tee -a ${logFile}
echo "Reference genome: " $alignRef | tee -a ${logFile}
echo "Barcode length: " $barcodeLength | tee -a ${logFile}
echo "Spacer length: " $spacerLength | tee -a ${logFile}
echo "Pre-tag_to_header read length: " $(($readLength+$barcodeLength+$spacerLength)) | tee -a ${logFile}
echo "Post-tag_to_header read length: " $readLength | tee -a ${logFile}
echo "Repetitive tag filter length: " $repFilt | tee -a ${logFile}
echo "Minimum family size: " $minMem | tee -a ${logFile}
echo "Maximum family size: " $maxMem | tee -a ${logFile}
echo "Consensus cutoff: " $cutOff | tee -a ${logFile}
echo "Consensus N cutoff: " $nCutOff | tee -a ${logFile}
echo "Read types: " $readTypes | tee -a ${logFile}
echo "Filters: " $filtersSet | tee -a ${logFile}
echo "Trim beginning: "$trimBeginningB"-"$trimBeginningE | tee -a ${logFile}
echo "Trim end: " $trimEndB"-"$trimEndE | tee -a ${logFile}
echo "Target region: " $TARGETREGION | tee -a ${logFile}
echo "SNPs to remove:" $SNPS2PURGE | tee -a ${logFile}
echo "" | tee -a ${logFile}


# Run tag_to_header.py on imput files
echo "Starting Run" | tee -a ${logFile}
echo "tag_to_header starting"  | tee -a ${logFile}
date | tee -a ${logFile}
echo "" | tee -a ${logFile}
python ${DSpath}/tag_to_header.py --infile1 $read1in --infile2 $read2in --outprefix ${runIdentifier} --tagstats


# Align sequences
echo "Aligning with BWA" | tee -a ${logFile}
date | tee -a ${logFile}
bwa aln $alignRef ${runIdentifier}.seq1.smi.fq > ${runIdentifier}.seq1.aln
bwa aln $alignRef ${runIdentifier}.seq2.smi.fq > ${runIdentifier}.seq2.aln
bwa sampe -s $alignRef ${runIdentifier}.seq1.aln ${runIdentifier}.seq2.aln ${runIdentifier}.seq1.smi.fq ${runIdentifier}.seq2.smi.fq > ${runIdentifier}.pe.sam


# Sort and index aligned sequences
echo "Sorting aligned raw sequences and indexing" | tee -a ${logFile}
date | tee -a ${logFile}
samtools view -Sbu ${runIdentifier}.pe.sam | samtools sort - ${runIdentifier}.pe.sort
samtools index ${runIdentifier}.pe.sort.bam


# Run Consensus Maker
echo "Starting Consensus Maker" | tee -a ${logFile}
date | tee -a ${logFile}
python ${DSpath}/ConsensusMaker.py --infile ${runIdentifier}.pe.sort.bam --tagfile ${runIdentifier}.pe.tagcounts --outfile ${runIdentifier}.sscs.bam --minmem $minMem --maxmem $maxMem --readlength $readLength --cutoff $cutOff --Ncutoff $nCutOff --read_type $readTypes --filt $filtersSet --isize $iSize


# Sort and index SSCSs
echo "Sorting SSCSs and indexing" | tee -a ${logFile}
date | tee -a ${logFile}
samtools view -bu ${runIdentifier}.sscs.bam | samtools sort - ${runIdentifier}.sscs.sort
samtools index ${runIdentifier}.sscs.sort.bam


# Run Duplex Maker
echo "Starting Duplex Maker" | tee -a ${logFile}
date  | tee -a ${logFile}
python ${DSpath}/DuplexMaker.py --infile ${runIdentifier}.sscs.sort.bam --outfile ${runIdentifier}.dcs.bam --Ncutoff $nCutOff --readlength $readLength


# Align DCSs
echo "Aligning DCSs" | tee -a ${logFile}
date | tee -a ${logFile}
bwa aln $alignRef ${runIdentifier}.dcs.r1.fq > ${runIdentifier}.dcs.r1.aln
bwa aln $alignRef ${runIdentifier}.dcs.r2.fq > ${runIdentifier}.dcs.r2.aln
bwa sampe -s $alignRef ${runIdentifier}.dcs.r1.aln ${runIdentifier}.dcs.r2.aln ${runIdentifier}.dcs.r1.fq ${runIdentifier}.dcs.r2.fq > ${runIdentifier}.dcs.aln.sam


# Sort and index aligned DCSs
echo "Sorting aligned DCSs and indexing" | tee -a ${logFile}
date | tee -a ${logFile}
samtools view -Sbu ${runIdentifier}.dcs.aln.sam | samtools sort - ${runIdentifier}.dcs.aln.sort
samtools index ${runIdentifier}.dcs.aln.sort.bam

———————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————

# Filter out unmapped reads & re-index
echo "Filtering unmapped DCSs and indexing" | tee -a ${logFile}
date | tee -a ${logFile}
samtools view -F4 ${runIdentifier}.dcs.aln.sort.bam | samtools view -Sb -T $alignRef - > ${runIdentifier}.dcs.filt.bam
samtools index ${runIdentifier}.dcs.filt.bam


# End-clip DCS reads
echo "End-clipping DCS reads and indexing" | tee -a ${logFile}
date | tee -a ${logFile}
java -jar -Xmx2g $PICARDPATH/AddOrReplaceReadGroups.jar INPUT=${runIdentifier}.dcs.filt.bam OUTPUT=${runIdentifier}.dcs.filt.readgroups.bam RGLB=Eeny RGPL=Meeny RGPU=Miny RGSM=Moe
samtools index ${runIdentifier}.dcs.filt.readgroups.bam
java -jar -Xmx8g $GATKPATH/GenomeAnalysisTK.jar -T ClipReads -I ${runIdentifier}.dcs.filt.readgroups.bam -o ${runIdentifier}.dcs.clipped.bam -R $alignRef --cyclesToTrim "$trimEndB"-"$trimEndE,$trimBeginningB"-"$trimBeginningE" --clipRepresentation SOFTCLIP_BASES --fix_misencoded_quality_scores
samtools index ${runIdentifier}.dcs.clipped.bam


# Locally realign clipped DCS reads (MUST be done for clipOverlap to work)
java -Xmx4g  -jar $GATKPATH/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $alignRef -I ${runIdentifier}.dcs.clipped.bam -o ${runIdentifier}.dcs.clipped.bam.readgroups.intervals
java -Xmx4g -jar $GATKPATH/GenomeAnalysisTK.jar -T IndelRealigner -R $alignRef -I ${runIdentifier}.dcs.clipped.bam -targetIntervals ${runIdentifier}.dcs.clipped.bam.readgroups.intervals -o ${runIdentifier}.dcs.clipped.realign.bam
samtools index ${runIdentifier}.dcs.clipped.realign.bam
rm ${runIdentifier}.dcs.clipped.bam.readgroups.intervals


# Clip overlapping nucleotides (softclipping, default seems to be clip from first read pair)
echo "Cliping overlapping nucleotides and recording stats" | tee -a ${logFile}
date | tee -a ${logFile}
bam clipOverlap --in ${runIdentifier}.dcs.clipped.realign.bam --out ${runIdentifier}.dcs.clipped.deoverlap.bam --stats 2> overlapping-reads.stats


# Filter out non-perfectly-mapping reads (those with mapping score less than 37), index. Note that perfect mapping in PE is maybe actually higher than 37? Need check. Different than in SE.
echo "Filtering out reads mapping with <q37 confidence and indexing" | tee -a ${logFile}
date | tee -a ${logFile}
samtools view -bh -q37 ${runIdentifier}.dcs.clipped.deoverlap.bam > ${runIdentifier}.dcs.clipped.deoverlap.q37.bam
samtools index ${runIdentifier}.dcs.clipped.deoverlap.q37.bam


# Create a more easily visualized bam file with only mutant reads included, index
echo "Creating mutation only BAM file" | tee -a ${logFile}
date | tee -a ${logFile}
samtools view ${runIdentifier}.dcs.clipped.deoverlap.q37.bam | awk '{if ($17 != "NM:i:0") print $0 }' - > ${runIdentifier}.dcs.clipped.deoverlap.q37.muts.sam
samtools view -bT $alignRef ${runIdentifier}.dcs.clipped.deoverlap.q37.muts.sam > ${runIdentifier}.dcs.clipped.deoverlap.q37.muts.bam
samtools index ${runIdentifier}.dcs.clipped.deoverlap.q37.muts.bam


# Create pileup. (Note: -Q0 filter is apparently critical for making mpileup work on softclipped reads)
echo "Creating pileup(s)" | tee -a ${logFile}
date | tee -a ${logFile}
samtools mpileup -Q0 -B -A -d 500000 -f $alignRef ${runIdentifier}.dcs.clipped.deoverlap.q37.bam > ${runIdentifier}.dcs.clipped.deoverlap.q37.pileup # (this is main one)
samtools mpileup -Q0 -B -A -d 500000 -f $alignRef ${runIdentifier}.dcs.clipped.bam > ${runIdentifier}.dcs.clipped.pileup # (this is from pre-overlap trim data to permit mutsbyposition plot) 


# Trim pileup to target region
echo "Trimming pileup to target region" | tee -a ${logFile}
date | tee -a ${logFile}
cat ${runIdentifier}.dcs.clipped.deoverlap.q37.pileup | bash $TARGETREGION > ${runIdentifier}.dcs.clipped.deoverlap.q37.region.pileup


# Generate pileup with SNP sites removed
echo "Removing SNPs from pileup" | tee -a ${logFile}
date | tee -a ${logFile}
cat ${runIdentifier}.dcs.clipped.deoverlap.q37.region.pileup | bash $SNPPURGE > ${runIdentifier}.dcs.clipped.deoverlap.q37.region.deSNP.pileup


# Generate final DCS mutation data files
echo "Generating DCS mutation files" | tee -a ${logFile}
date | tee -a ${logFile}
# Below depth -d is 1 (count all sites regardless of depth) and -C is 1 (turn off clonality cutoff). Default N filter is 5%. -u is counting uniques only .
cat ${runIdentifier}.dcs.clipped.deoverlap.q37.region.deSNP.pileup | python $DSpath/CountMuts.py -d -C 1 -u > ${runIdentifier}.dcs.clipped.deoverlap.q37.region.deSNP.pileup.countmuts
# Below depth -d is 1 (count all sites regardless of depth) and -C is 1 (turn off clonality cutoff). N filter 1 (no max). No -u (keep all, not just uniques) 
cat ${runIdentifier}.dcs.clipped.deoverlap.q37.region.deSNP.pileup | python $DSpath/mut-position.py -d 1 -n 1 > ${runIdentifier}.DCS-muts.noSNP.txt
cat ${runIdentifier}.dcs.clipped.deoverlap.q37.region.deSNP.pileup | python $DSpath/mut-position.py > ${runIdentifier}.dcs.clipped.deoverlap.q37.region.deSNP.pileup.mutpos


# Generate pileup with only SNP sites of interest and create mutation file from this
echo "Generating DCS SNP files" | tee -a ${logFile}
date | tee -a ${logFile}
FILE2SIFT=./${runIdentifier}.dcs.clipped.deoverlap.q37.region.pileup; export FILE2SIFT	         ### This variable must be defined for called script to work
bash $SNPSIFT > ${runIdentifier}.dcs.clipped.deoverlap.q37.region.SNP-only.pileup
cat ${runIdentifier}.dcs.clipped.deoverlap.q37.region.SNP-only.pileup | python $DCSPATH/mut-position.py -d 1 -n 1 > ${runIdentifier}.DCS-SNP-only.txt


# Generate statistics files:
echo "Generating data stats" | tee -a ${logFile}
date | tee -a ${logFile}

python $DSpath/muts_by_read_position.py --infile ${runIdentifier}.dcs.clipped.pileup --outfile ${runIdentifier}.dcs.clipped.mutsbyreadposition.png --rlength 85
	# (Note this is from DCS data after clipping but before any overlap trim, q37 filtering, target region or SNP filtering)

bash $DSpath/get-flagStats-PE.txt > flagstats.stats
	# (Number of reads at different steps)

bash $DSpath/get-idxStats-PE.txt > idxStats.stats
	# (Distribution chromosome mapping at different steps)

samtools view ${runIdentifier}.sscs.sort.bam | cut -f9 | awk '{count[$1]++} END {for(j in count) print j," \t", count[j]}' | sort -n | awk '{if (($1 > 0) && ($1 < 800)) print $1 "\t" $2}' > iStats.pe.txt
	# (iStats for raw reads)

samtools view ${runIdentifier}.dcs.clipped.deoverlap.q37.bam | cut -f9 | awk '{count[$1]++} END {for(j in count) print j," \t", count[j]}' | sort -n | awk '{if (($1 > 0) && ($1 < 800)) print $1 "\t" $2}' > iStats.dcs.final.txt
	# (iStats for final DCS reads)

cat ${runIdentifier}.dcs.clipped.deoverlap.q37.region.deSNP.pileup | cut -f4 | awk '{count[$1]++} END {for(j in count) print j," \t",count[j]}' | sort -nr > DCS-final-depth-stats.txt
	# (Distribution of duplex depth. First column is depth at position, 2nd is # of positions with that depth)


# Cleanup workspace and finish
rm *.sam
echo "Finished with run " $runIdentifier | tee -a ${logFile}
date | tee -a ${logFile}









