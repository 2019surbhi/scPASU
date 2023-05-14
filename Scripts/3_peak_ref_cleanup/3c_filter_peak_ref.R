#!/usr/bin/env Rscript

library(dplyr)
library(data.table)

source('/home/sonas/beegfs/APA/scPASU/scripts/scPASU_functions.R')

args <- commandArgs(trailing = TRUE)


ref_file=args[1]
cov_mat_file=args[2]
out=args[3]
frefix=args[4]
min_cov=as.numeric(args[5])
is_minus=as.logical(args[6]) # TRUE or FALSE


# HPC paths

ref_file<-'/home/sonas/beegfs/APA/scPASU/output/3_RefinePeakRef/3a_assign_TU/u10_uro_peak_universe_updated.txt'
cov_mat_file<-'/home/sonas/beegfs/APA/scPASU/output/3_RefinePeakRef/3b_PeakCoverage2/u10_uro_minus_filtered_peak_count.rds'

# Read inputs
ref<-fread(ref_file)
cov_mat<-readRDS(cov_mat_file)

cov_mat_fil<-cov_mat[which(cov_mat$pcov_pct>=min_cov),]

cat('dims: \n')
dim(cov_mat)
dim(cov_mat_fil)

keep<-match(cov_mat_fil$anno,ref$final_annotation)
ref_filtered<-ref[keep,]

# Subset TUs that have now single peaks
cnt<-ref_filtered %>% group_by(tu) %>% tally()
cnt<-cnt %>% arrange(.,n)
r<-which(cnt$n==1)
p0<-ref_filtered[ref_filtered$tu %in% cnt$tu[r],]
pn<-ref_filtered[!(ref_filtered$tu %in% cnt$tu[r]),]

pn_updated<-create_final_annotation_col(pn,is_minus=is_minus)
ref_filtered_updated<-rbind(p0,pn_updated)


write.table(ref_filtered_updated,paste0(out,frefix,'_filtered_peak_ref.txt'),col.names=TRUE,row.names=FALSE,sep='\t',quote=FALSE)
