#!/usr/bin/env Rscript

args<-commandArgs(trailing=TRUE)
stats_dir=args[1]
outdir=args[2]
fprefix=args[3]
strand=args[4]

library(data.table)
library(dplyr)

cat('stats dir: ', stats_dir, '\n')
cat('strand: ', strand, '\n')

setwd(stats_dir)

# read files
peaks_file<-fread(paste0(stats_dir,'total_peaks_stats_',strand,'.txt'))
polya_peaks_file<-fread(paste0(stats_dir,'peaks_supported_by_polya_',strand,'.txt'))
non_polya_peaks_file<-fread(paste0(stats_dir,'peaks_not_supported_by_polya_',strand,'.txt'))
min_polya_peaks_file<-fread(paste0(stats_dir,'peaks_supported_by_min_polya_',strand,'.txt'))

polya_reads_file<-fread(paste0(stats_dir,'total_polya_stats_',strand,'.txt'))
polya_peaks_reads_file<-fread(paste0(stats_dir,'polya_reads_supporting_peaks_',strand,'.txt'))
min_polya_peaks_reads_file<-fread(paste0(stats_dir,'min_polya_reads_supporting_peaks_',strand,'.txt'))

# Extract stats and merge in a single df

nms<-peaks_file[,2] %>% pull()
n<-length(nms)

#Remove 'total'
nms<-nms[-n]

# String split by _ and then extract 2nd sub element
nms<-strsplit(nms,split='_')
chr<-sapply(nms,'[[',2)

chr<-c(chr,'total')
peaks<-peaks_file[,1]
polya_peaks<-polya_peaks_file[,1]
non_polya_peaks<-non_polya_peaks_file[,1] %>% pull()
min_polya_peaks<-min_polya_peaks_file[,1]

polya_reads<-polya_reads_file[,1]
polya_peaks_reads<-polya_peaks_reads_file[,1] %>% pull() #
min_polya_peaks_reads<-min_polya_peaks_reads_file[,1] %>% pull()


# Find sum for the tables that don't have sum calculated
total<-sum(non_polya_peaks)
non_polya_peaks<-c(non_polya_peaks,total)

total<-sum(polya_peaks_reads)
polya_peaks_reads<-c(polya_peaks_reads,total)


# Merge
stats<-cbind(peaks,polya_peaks,non_polya_peaks,min_polya_peaks,polya_reads,polya_peaks_reads,min_polya_peaks_reads)
colnames(stats)<-c('total_peaks','polya_peaks','non_polya_peaks','min_polya_peaks','polya_reads','polya_reads_supporting_peaks','min_polya_reads_supporting_peaks')
stats<-cbind(chr,stats)

# Save
write.csv(stats,paste0(outdir,fprefix,'_',strand,'_stats.csv'))
