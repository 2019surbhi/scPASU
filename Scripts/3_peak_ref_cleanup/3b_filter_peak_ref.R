#!/usr/bin/env Rscript

library(dplyr)
library(data.table)

args <- commandArgs(trailing = TRUE)


pcov_file=args[1]
out=args[2]
fprefix=args[3]
min_cov=as.numeric(args[4])

# HPC paths
pcov_file<-'/home/sonas/beegfs/APA/scPASU/output/3_RefinePeakRef/3b_PeakCoverage2/u10_uro_minus_peak_count_updated.rds'
out<-'/home/sonas/beegfs/APA/scPASU/output/3_RefinePeakRef/3c_updated_peaks/'
fprefix<-'u10_uro_minus'

pcov_mat<-readRDS(pcov_file)
rem<-which(pcov_mat$pcov_pct<min_cov)
pcov_mat2<-pcov_mat[-rem,]
