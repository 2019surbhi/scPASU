#!/usr/bin/env Rscript

library(dplyr)
library(data.table)

args <- commandArgs(trailing = TRUE)

ref_file=args[1]
counts_file=args[2]

# HPC paths 
#ref_file='/home/sonas/beegfs/APA/scPASU/output/3_RefinePeakRef/3a_assign_TU/u10_uro_peak_universe_minus_updated.saf'
#counts_file='/home/sonas/beegfs/APA/scPASU/output/3_RefinePeakRef/3b_PeakCoverage2/u10_uro_minus_peak_by_cell_count.rds'


#Read inputs
ref<-fread(ref_file)
counts<-readRDS(counts_file)
counts_mat<-counts$counts


tu<-strsplit(ref$GeneID,split=':') %>% sapply(.,'[[',1)
gene<-strsplit(ref$GeneID,split=':') %>% sapply(.,'[[',2)
peak<-strsplit(ref$GeneID,split=':') %>% sapply(.,'[[',3)

ref$tu<-tu
ref$gene<-gene
ref$peak<-peak

updated_mat<-data.frame()

for(i in 1:length(tu))
{
r<-grep(tu[i],rownames(counts_mat))
sub<-counts_mat[r,] %>% as.data.frame()
cov<-sum(sub)
pcov<-apply(sub,1,sum)
pcov_pct<-(pcov/cov)*100
sub<-cbind(sub,pcov_pct)

if(nrow(updated_mat)==0)
{
  updated_mat<-sub
}else{
  updated_mat<-rbind(updated_mat,sub)
}
}

colnames(updated_mat)<-c('counts','cov_pct')
out<-dir(counts_file)
fname<-gsub('.rds','_updated.rds',basename(counts_file))
saveRDS(updated_mat,paste0(out,fname))
