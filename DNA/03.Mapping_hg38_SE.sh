#!/bin/bash

outdir=$1
fq=$2
bismark --bowtie2 -p 24 -quiet --non_directional --unmapped --temp_dir $outdir -o $outdir /path/to/WholeGenomeFasta $fq
