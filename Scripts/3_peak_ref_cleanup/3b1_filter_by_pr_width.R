#!/usr/bin/env Rscript

library(dplyr)
library(data.table)

args <- commandArgs(trailing = TRUE)


ref_file=args[1]
out=args[3]
frefix=args[4]
max_width=as.numeric(args[5])

# HPC paths

#ref_file<-'/home/sonas/beegfs/APA/scPASU/output/3_RefinePeakRef/3a_assign_TU/u10_uro_peak_universe_updated.txt'

# Read inputs
ref<-fread(ref_file)
ref_filtered<-ref[which(ref$pr_width<=max_width),]

write.table(ref_filtered,paste0(out,frefix,'_filtered_peak_ref.txt'),col.names=TRUE,row.names=FALSE,sep='\t',quote=FALSE)

