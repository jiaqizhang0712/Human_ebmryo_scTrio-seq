#!/bin/bash

outdir=$1
fq1=$2
fq2=$3

bismark --bowtie2 -p 24 -quiet --non_directional --unmapped --temp_dir $outdir -o $outdir /path/to/WholeGenomeFasta  -1 $fq1 -2 $fq2