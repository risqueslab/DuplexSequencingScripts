#!/bin/bash

#  Unified_pipeline.sh
#  
#
#  Created by Kennedy Lab on 5/17/16.
#
alignRef='/Users/kennedylab/Bioinformatics/genomes/Human_Genome_Sequence/human_g1k_v37.fasta'


folderList='CODIS1 CODIS2 CODIS3 P53_sonicated_1 P53_sonicated_2 P53_sonicated_3 P53_cas9_1 P53_cas9_2 P53_cas9_3'
for elmt in $folderList;
do
    cd ${elmt}
    java -jar ~/Bioinformatics/programs/picard-tools-2.2.4/picard.jar FastqToSam F1=${elmt}.seq1.fastq.gz F2=${elmt}.seq2.fastq.gz O=/dev/stdout SM=${elmt}|python ~/Bioinformatics/programs/Duplex-Sequencing-master/UnifiedConsensusMaker.py --input /dev/stdin --taglen 10 --spacerlen 1 --write-sscs --prefix ${elmt} --tagstats

    bwa aln $alignRef ${elmt}_read1_sscs.fq.gz > ${elmt}.sscs.read1.aln
    bwa aln $alignRef ${elmt}_read2_sscs.fq.gz > ${elmt}.sscs.read2.aln
    bwa sampe -s $alignRef ${elmt}.sscs.read1.aln ${elmt}.sscs.read2.aln ${elmt}_read1_sscs.fq.gz ${elmt}_read2_sscs.fq.gz | samtools sort -o ${elmt}.sscs.sort.bam -

    bwa mem $alignRef ${elmt}_read1_sscs.fq.gz ${elmt}_read2_sscs.fq.gz|samtools sort -o ${elmt}_mem.sscs.sort.bam -

    bwa aln $alignRef ${elmt}_read1_dcs.fq.gz > ${elmt}.dcs.read1.aln
    bwa aln $alignRef ${elmt}_read2_dcs.fq.gz > ${elmt}.dcs.read2.aln
    bwa sampe -s $alignRef ${elmt}.dcs.read1.aln ${elmt}.dcs.read2.aln ${elmt}_read1_dcs.fq.gz ${elmt}_read2_dcs.fq.gz | samtools sort -o ${elmt}.dcs.sort.bam -
    bwa mem $alignRef ${elmt}_read1_dcs.fq.gz ${elmt}_read2_dcs.fq.gz|samtools sort -o ${elmt}_mem.dcs.sort.bam -

    samtools index ${elmt}.sscs.sort.bam
    samtools index ${elmt}_mem.sscs.sort.bam
    samtools index ${elmt}.dcs.sort.bam
    samtools index ${elmt}_mem.dcs.sort.bam
    bash ~/Bioinformatics/programs/Duplex-Sequencing-master/PostDCSProcessing_human.sh ${elmt}_mem.dcs.sort.bam ~/Bioinformatics/genomes/Human_Genome_Sequence/human_g1k_v37.fasta 100 0.00 0.01
    cd ..
done
