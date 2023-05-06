#!/usr/bin/env Rscript

library(data.table)
library(dplyr)


file=args[1]
out=args[2]
fname=args[3]

ref<-fread(file)
rem<-which(ref$V5!=ref$V11)
ref2<-ref[-rem,]
write.table(ref2,paste0(out,fname),quote=FALSE,row.names=FALSE,colnames=FALSE,sep='\t')

