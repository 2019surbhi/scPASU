#!/usr/bin/env Rscript

library(dplyr)
library(data.table)

args <- commandArgs(trailing = TRUE)
bed_file=args[1]
outpath=args[2]
min_polya=as.numeric(args[3])

cat('Reading bed file: ',bed_file, '\n')
tab<-fread(bed_file)

cat('dim of original: ', dim(tab),'\n')

keep<-which(tab$polya_count>=min_polya)
tab2<-tab[keep,]

cat('dim of filtered: ', dim(tab2),'\n')

fname<-gsub('.bed',paste0('_',min_polya,'polya.bed'),basename(bed_file))

cat('Writing file: ',fname,' to path: ',outpath,'\n')
write.table(tab2,paste0(outpath,fname),sep='\t')

