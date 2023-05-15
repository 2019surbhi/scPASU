#!/usr/bin/env Rscript

library(dplyr)
library(data.table)

args <- commandArgs(trailing = TRUE)


ref_file=args[1]
cov_mat_file[2]
out=args[3]
frefix=args[4]
mincov=as.numeric(args[5])

ref<-fread(ref_file)
cov_mat<-readRDS(cov_mat_file)
sub<-cov_mat[which(cov_mat$pcov_pct>=mincov),]

select<-match(sub$peakID,ref$peakID)
ref_filtered<-ref[select,]

write.table(ref_filtered,paste0(out,frefix,'_cov_filtered_peak_ref.txt'),col.names=TRUE,row.names=FALSE,sep='\t',quote=FALSE)

