#!/usr/bin/env Rscript

library(data.table)
library(dplyr)

args<-commandArgs(trailingOnly = TRUE)

ref_file_path=args[1]
out=args[2]
prefix=args[3]

# Read ref file #
ref_file<-read.table(ref_file_path,header=TRUE)

# Select relevant columns
cols<-c('final_annotation','chr','start','end','strand')
select<-match(cols,colnames(ref_file))
saf_ref<-ref_file[,select]
colnames(saf_ref)<-c('GeneID','Chr','Start','End','Strand') 

cat('Creating SAF ref file \n')
write.table(saf_ref,paste0(out,prefix,'_updated.saf'),sep='\t',quote=FALSE,row.names=FALSE)
