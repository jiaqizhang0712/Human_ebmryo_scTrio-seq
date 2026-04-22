args<-commandArgs(TRUE)
suppressMessages(library(HMMcopy))

sample.name <- args[1] # sample name
rfile <- args[2] # .wig file 
window <- "1000000"

gfile <- paste0("/path/to/HMMcopy/hmmcopy_utils-master/hg38/hg38.gc.",window,".wig") # gc file produced by hmmcopy toolkit
mfile <- paste0("/path/to/HMMcopy/hmmcopy_utils-master/hg38/hg38.map.",window,".bw") # map file produced by hmmcopy toolkit

normal_reads <- wigsToRangedData(rfile, gfile, mfile)
normal_copy <- correctReadcount(normal_reads)
normal_copy$chr <- as.factor(normal_copy$chr)
tumour_segments <- HMMsegment(normal_copy)

write.csv(normal_copy,file = paste0(sample.name,".hmmcopy.win",window,".corrected.csv"),row.names = F,quote = F)
write.csv(tumour_segments$segs,file = paste0(sample.name,".hmmcopy.win",window,".segment.csv"),row.names = F,quote = F)