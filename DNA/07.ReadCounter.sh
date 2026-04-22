outdir=$1 # output dir
sample=$2 # sample name
bam=$3 # rmdup bam file 
window=$4 # window size, defualt is 1000000

readCounter -w $window -c chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY $bam > ${outdir}/${sample}.window${window}.readcounts.wig