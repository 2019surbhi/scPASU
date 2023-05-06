#!/usr/bin/env Rscript

library(data.table)
library(readr)
library(dplyr)

args <- commandArgs(trailing = TRUE)
dir=args[1]
bc_dir=args[2]
sample=args[3]

#dir<-'/home/sonas/beegfs/APA/scAPA/ureter10/peak_cell_mat_test/'

#bc_dir<-paste0(dir,'uro_barcodes/')
bc_file<-fread(paste0(bc_dir,sample,'_uro_barcodes.tsv'), header=FALSE)
bc<-bc_file$V1


out<-paste0(dir,'/split_barcodes/',sample,'/')

#for(i in 1:length(bc))
#{
# write_tsv(bc[i],paste0(out,bc[i],'.tsv'),col_names=FALSE)
#}


lapply(1:length(bc), function(x){write_tsv(as.data.frame(bc[x]),paste0(out,bc[x],'.tsv'),col_names=FALSE,quote="none")})

