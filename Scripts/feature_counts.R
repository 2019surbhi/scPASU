#!/usr/bin/env Rscript

library(data.table)
library(dplyr)
library(Rsubread)


args<-commandArgs(trailingOnly = TRUE)


bam_file<-args[1]
ref<-args[2] # GTF of saf format features ref
out<-args[3]
prefix<-args[4]
cores<-args[5]
isGTF<-args[6] # Enter yes or no


if(isGTF=='yes')
{
  counts<-featureCounts(bam_file, annot.ext=ref, isGTFAnnotationFile=TRUE, GTF.attrType='gene_name', strandSpecific=1, largestOverlap = T, nthreads = cores)
  fname<-paste0(out,prefix,'_gene_count.rds')
}else{
  counts<-featureCounts(bam_file, annot.ext = ref, isGTFAnnotationFile = FALSE,strandSpecific=1, largestOverlap = T, nthreads = cores)
  fname<-paste0(out,prefix,'_peak_count.rds')
}

saveRDS(counts, fname)
