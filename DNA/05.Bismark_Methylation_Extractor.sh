#!/bin/bash

genome_folder=/path/to/WholeGenomeFasta
bam=$1
outdir=$2

######  For CpG only  08.3.bismark_methylation_extractor.CpG ######
bismark_methylation_extractor  --multicore 8 --gzip --cytosine_report --report --comprehensive --bedGraph --genome_folder $genome_folder -o $outdir $bam