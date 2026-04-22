#!/bin/bash

PICARD="java -Xmx4g -jar /share/home/zhangjiaqi/software/picard-tools-1.141/picard.jar"

sampledir=$1
sample=$2

# Input BAM files:
#   bam1: paired-end Bismark alignment BAM
#   bam2: single-end Bismark alignment BAM for read1 unmapped reads
#   bam3: single-end Bismark alignment BAM for read2 unmapped reads
bam1=$3 
bam2=$4
bam3=$5


### Sort input BAMs
samtools sort -o ${bam1}.sort.bam $bam1
samtools sort -o ${bam2}.sort.bam $bam2
samtools sort -o ${bam3}.sort.bam $bam3

### Index sorted BAMs
samtools index ${bam1}.sort.bam
samtools index ${bam2}.sort.bam
samtools index ${bam3}.sort.bam

### Merge BAMs
samtools merge ${sampledir}/${sample}.merged.bam  ${bam1}.sort.bam  ${bam2}.sort.bam  ${bam3}.sort.bam
samtools index ${sampledir}/${sample}.merged.bam

### Mark and remove duplicates 
$PICARD MarkDuplicates I=${sampledir}/${sample}.merged.bam O=${sampledir}/${sample}.merged.rmdup.bam M=${sampledir}/${sample}.merged.rmdup.txt  REMOVE_DUPLICATES=true 2> ${sampledir}/${sample}.picard.log

samtools index ${sampledir}/${sample}.merged.rmdup.bam


### Filter reads by mapping quality, q10
samtools view -h -q 10 -o ${sampledir}/${sample}.merged.rmdup.q10.bam  ${sampledir}/${sample}.merged.rmdup.bam 2>${sampledir}/${sample}.q10.log

samtools index ${sampledir}/${sample}.merged.rmdup.q10.bam 