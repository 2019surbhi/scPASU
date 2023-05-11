#!/usr/bin/env Rscript

library(dplyr)
library(data.table)

args <- commandArgs(trailing = TRUE)


pcov_file=args[1]
out=args[2]
fprefix=args[3]

pcov_mat<-readRDS(pcov_file)


# HPC paths
fprefix<-'u10_uro_minus'
out<-'/home/sonas/beegfs/APA/scPASU/output/3_RefinePeakRef/3c_updated_peaks/
pcov_file<-'/home/sonas/beegfs/APA/scPASU/output/3_RefinePeakRef/3b_PeakCoverage2/u10_uro_minus_peak_count_updated.rds'
