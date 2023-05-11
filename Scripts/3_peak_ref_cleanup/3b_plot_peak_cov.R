#!/usr/bin/env Rscript

library(dplyr)
library(data.table)

source('/home/sonas/beegfs/APA/scPASU/scripts/scPASU_functions.R')

args <- commandArgs(trailing = TRUE)

pcov_file=args[1]
out=args[2]
fprefix=args[3]
w=args[4]
h=args[5]
x=args[6]
y=args[7]

#HPC paths
# pcov_file='u10_uro_plus_peak_by_cell_count_updated.rds'
# out='/home/sonas/beegfs/APA/scPASU/output/3_RefinePeakRef/3b_PeakCoverage2/'
# fprefix='u10_uro_plus'

pcov_mat<-readRDS(pcov_file)

dat<-pcov_mat$pcov_pct

options(bitmapType='cairo')
pdf(paste0(out,fprefix,'_peak_cov_pct.pdf'),width=w, height=h)
get_histogram(dat=dat,binw=1,col='royal blue',fill='salmon', x_lim=c(0,100),y_lim=c(0,5000),title='Peak read count percentage',x_lab='peak_read_pct')
dev.off()

