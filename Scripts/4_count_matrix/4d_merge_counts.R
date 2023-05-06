#!/usr/bin/env Rscript

library(data.table)
library(readr)
library(dplyr)

args<-commandArgs(trailingOnly = TRUE)
counts_path<-args[1]
file_prefix<-args[2]
out_dir<-args[3]
  
setwd(counts_path)

f<-list.files(counts_path,pattern = '.rds',full.names = TRUE)
dat<-lapply(f,readRDS)

# The outputs of SubRead are stored as list so dat is list of lists - extracting features from one of the countr matrix
counts<-rownames(dat[[1]]$counts)

 for(i in 1:length(dat))
 {
  cnt<-dat[[i]]$counts
  counts<-cbind(counts,cnt)
  
  }

colnames(counts)<-gsub('.bam','',colnames(counts))

counts<-as.data.frame(counts)
counts<-counts[,-1]

write.table(counts,paste0(out_dir,file_prefix,'_counts.txt'),row.names=TRUE,col.names=TRUE,sep='\t')
 
 
