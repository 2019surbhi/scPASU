#!/usr/bin/env Rscript

library(dplyr)
library(data.table)

source('/home/sonas/beegfs/APA/scPASU/scripts/scPASU_functions.R')

args <- commandArgs(trailing = TRUE)


pcov_file=args[1]
out=args[2]
fprefix=args[3]
min_cov=as.numeric(args[4])
#is_minus=as.logical(args[5]) # TRUE or FALSE

ref<-fread(ref_file)

keep<-match(cov_mat_fil$anno,ref$final_annotation)
ref_filtered<-ref[keep,]

# Subset TUs that have now single peaks
#cnt<-ref_filtered %>% group_by(tu) %>% tally()
#cnt<-cnt %>% arrange(.,n)
#r<-which(cnt$n==1)
#p0<-ref_filtered[ref_filtered$tu %in% cnt$tu[r],]
#pn<-ref_filtered[!(ref_filtered$tu %in% cnt$tu[r]),]

#pn_updated<-create_final_annotation_col(pn,is_minus=is_minus)
#ref_filtered_updated<-rbind(p0,pn_updated)

#write.table(ref_filtered_updated,paste0(out,frefix,'_filtered_peak_ref.txt'),col.names=TRUE,row.names=FALSE,sep='\t',quote=FALSE)

write.table(ref_filtered,paste0(out,frefix,'_filtered_peak_ref.txt'),col.names=TRUE,row.names=FALSE,sep='\t',quote=FALSE)


