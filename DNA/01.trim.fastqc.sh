#!/bin/bash

outdir=$1
fq1=$2
fq2=$3

# TrimGalore-0.6.10
trim_galore --fastqc --clip_R1 9 --clip_R2 9  --quality 20 --phred33 --stringency 3 --gzip --length 30 --paired --retain_unpaired  --output_dir $outdir  $fq1 $fq2