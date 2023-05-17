#!/usr/bin/env Rscript

library(data.table)
library(dplyr)

argv<-commandArgs(trailing = TRUE)

counts_dir=argv[1]
fprefix=argv[2] 
out=argv[3]


f<-list.files(counts_dir,full.names = TRUE)
peak_counts<-lapply(f,read.table,header=TRUE,sep='\t')


# Replace . with - in count mat 
for(i in 1:length(peak_counts))
{
  colnames(peak_counts[[i]])<-gsub('.','-',colnames(peak_counts[[i]]),fixed = TRUE)
  colnames(peak_counts[[i]])<-gsub('uro_','',colnames(peak_counts[[i]]))
}

# Merge counts
merged_counts<-do.call(cbind,peak_counts)

write.table(merged_counts,paste0(out,fprefix,'_counts.txt'),sep='\t')
