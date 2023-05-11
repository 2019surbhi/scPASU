#!/usr/bin/env Rscript

library(dplyr)
library(data.table)

source('/home/sonas/beegfs/APA/scPASU/scripts/scPASU_functions.R')

args <- commandArgs(trailing = TRUE)


pcov_file=args[1]
out=args[2]
fprefix=args[3]
min_cov=as.numeric(args[4])
is_minus=as.logical(args[5]) # TRUE or FALSE

pcov_mat<-readRDS(pcov_file)
rem<-which(pcov_mat$pcov_pct<min_cov)
pcov_mat2<-pcov_mat[-rem,]

# Subset TUs that have now single peaks
cnt<-pcov_mat2 %>% group_by(tu) %>% tally()
cnt<-cnt %>% arrange(.,n)
r<-which(cnt$n==1)
p0_mat<-pcov_mat2[pcov_mat2$tu %in% cnt$tu[r],]
pn_mat<-pcov_mat2[!(pcov_mat2$tu %in% cnt$tu[r]),]

pn_mat2<-create_final_annotation_col(pn_mat,is_minus=is_minus)
pcov_mat3<-rbind(p0_mat,pn_mat)
saveRDS(pcov_mat3,paste0(out,fprefix,'_filtered_peak_count.rds'))

