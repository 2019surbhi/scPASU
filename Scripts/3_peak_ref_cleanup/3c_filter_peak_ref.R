#!/usr/bin/env Rscript

library(dplyr)
library(data.table)

source('/home/sonas/beegfs/APA/scPASU/scripts/scPASU_functions.R')

args <- commandArgs(trailing = TRUE)


ref_file=args[1]
cov_mat_file=args[2]


# HPC paths

ref_file<-'/home/sonas/beegfs/APA/scPASU/output/3_RefinePeakRef/3a_assign_TU/u10_uro_peak_universe_updated.txt'
cov_mat_file<-'/home/sonas/beegfs/APA/scPASU/output/3_RefinePeakRef/3b_PeakCoverage2/u10_uro_minus_filtered_peak_count.rds'

# Read inputs
ref<-fread(ref_file)
cov_mat<-readRDS(cov_mat_file)

