#!/usr/bin/env Rscript

library(dplyr)
library(data.table)

args <- commandArgs(trailing = TRUE)


ref_file=args[1]
out=args[2]
fprefix=args[3]


ref<-fread(ref_file)
ref<-as.data.frame(ref)

cols<-c('chr','start','end','final_annotation')
bed<-ref[,cols]
bed$score<-rep(0,nrow(bed))
bed$strand<-ref$strand

write.table(bed, paste0(out,fprefix,'.bed'),sep='\t', quote=FALSE, col.names=FALSE, row.names=FALSE)




