#!/usr/bin/env Rscript

library(dplyr)
library(stringi)
library(stringr)
library(readxl)
library(openxlsx)
library(data.table)
library(matrixStats)
library(gtools)



library(ggplot2)
library(reshape)

library(Rsamtools)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(edgeR)
library(DEXSeq)
library(rtracklayer)
library(corrplot)

source('scPASU_APA_testing_functions.R')

# User inputs

argv<-commandArgs(trailing = TRUE)

inputdir=argv[1] # full path to dir where merged count matrix is saved
outdir=argv[2] # path to output dir
counts_file=argv[3] # count matrix file name
meta_file=argv[4] # name of metadata file (saved from Seurat object, or a table which has cell cluster identity assigned to each cell)
peak_ref_file=argv[5] # peak file name with full path

# paths on HPC #

#inputdir<-'/home/sonas/beegfs/APA/scPASU/output/5_APA_testing/inputs/'
#outdir<-'/home/sonas/beegfs/APA/scPASU/output/5_APA_testing/outputs/'
#counts_file<-'u10_uro_counts.txt'
#meta_file<-'2021_04_01_ureter10_uro_PC50_res0.2_meta.xlsx'
#peak_ref_file<-'/home/sonas/beegfs/APA/scPASU/output/5_APA_testing/inputs/u10_uro_APA_testing_peak_ref.txt'

# Read merged counts
merged_counts<-read.table(paste0(inputdir,counts_file),sep='\t')
colnames(merged_counts)<-gsub('.','-',colnames(merged_counts),fixed=TRUE)

# Read meta data [save meta data from Seurat object - ensure the set of cells match that of the count matrix]
meta<-read.xlsx(paste0(inputdir,'2021_04_01_ureter10_uro_PC50_res0.2_meta.xlsx'))

### Load previous results ###

## Peak ref is jtu$join and final_annotation is the peak column 
peak_ref<-read.table(peak_ref_file,header=TRUE,sep='\t')

# Replace peak name column with final annotation as this is more meaningful name 
col<-which(colnames(peak_ref)=='final_annotation')
peak_ref$peak<-peak_ref[,col]
peak_ref<-peak_ref[,-col]

## Remove P0 peaks ##
plist<-strsplit(peak_ref$peak,split=':') %>% sapply(.,'[[',3)
rem<-which(plist=='P0') # 8066 P0 peaks
peak_ref<-peak_ref[-rem,]
 
# Read sample and comparison tables #

comp_file<-paste0(inputdir,'comp_group.csv')
compdt <- fread(comp_file,header=TRUE)
comps <- data.table(a=compdt$group1,b=compdt$group2)

# if the run fails run 5 or 10 comparisons at at a time
#bc<-comps
#comps<-comps[1:5,]

compnames <- paste0(comps$a,"_vs_",comps$b)

# Create stab

group<-c(comps$a,comps$b) %>% unique() %>% sort()
sample<-c(paste0(group,'_1'),paste0(group,'_2')) 
group<-rep(group,2)
n<-length(group)/2
stab<-cbind(sample,group) %>% as.data.frame()

stab$group<- stab$group %>% as.factor()

cov<-10
a<-comps$a
b<-comps$b

min_peak<-1
#adjust.var<-'batch'
ncpu<-20

# Per cluster input
peak_mat_lst<-list()
meanmat_lst<-list()

for(i in 1:nrow(comps))
{
c1<-as.numeric(gsub(".*?([0-9]+).*", "\\1", a[i]))
c2<-as.numeric(gsub(".*?([0-9]+).*", "\\1", b[i]))
clus_lst<-create_apa_inputs(merged_counts,meta,c1=c1,c2=c2,c1_nm=a[i],c2_nm=b[i],out=outdir,plot_corr = TRUE)

peak_mat_lst[[i]]<-clus_lst[[1]]
meanmat_lst[[i]]<-clus_lst[[2]]

}

# Order stab
#stab<-stab[match(colnames(peak_mat),stab$sample),]


# 1. DEXseq APA testing #

dotests <- function(a,b,ncpu=4,min_peak=1,meanmat,peak_mat)
{
  message("Testing: ",a," vs ",b)
  mycomp <- paste0(a,"_vs_",b)
  # Run without batch correction
  ta <- testAPA_sc(stab=stab,meanmat=meanmat,peak_ref=peak_ref,peak_mat=peak_mat,a=a,b=b,adjust.var=NULL,ncpu=ncpu,min_peak=min_peak)
  
  # Run with batch correction
  # taadj <- testAPA_sc(stab=stab,meanmat=meanmat,peak_ref=peak_ref,peak_mat=peak_mat,a=a,b=b,adjust.var=adjust.var,ncpu=ncpu,min_peak=min_peak)
  # list(ta=ta,taadj=taadj)
  
}

compnames <- paste0(comps$a,"_vs_",comps$b)
apa <- lapply(1:nrow(comps),function(x) dotests(a=comps[x,]$a,b=comps[x,]$b,min_peak=min_peak,meanmat = meanmat_lst[[x]],peak_mat = peak_mat_lst[[x]]))

names(apa) <- compnames


# Put the adj info into the main res so get one ta object for each

# tas <- lapply(apa,function(x){
#   x$res$noadj_p <- x$ta$res$p
#   x$res$noadj_padj <- x$ta$res$padj
#   
#   # x$taadj$res$batchadj_p <- x$taadj$res$p
#   # x$taadj$res$batchadj_padj <- x$taadj$res$padj
#    
#   x$res[,p:=NULL]
#   x$res[,padj:=NULL]
#   
#   return(x)
# })

#apa<-taadj

#apa <- tas

#output_file<-paste0(outdir,'apa.rd')
#save(apa,compress=T,file=output_file)


## 2. Run t test ##
nrepa<-2
nrepb<-2

apa.tt <- lapply(apa,testTTest,nrepa,nrepb)

#write.xlsx(apa.tt$res,paste0(outdir,'Differential_usage_ttest.xlsx'))

## 3. Run EdgeR ##
edger_prefix <- file.path(outdir,'edger_')
apa.alt <- lapply(1:length(apa.tt),function(x){testEdgeR(apa.tt[[x]],a=comps$a[x],comps$b[x])})

#write.xlsx(apa.alt$res,paste0(outdir,'Differential_usage_edger.xlsx'))
#save(apa.alt,file=paste0(outdir,'apa.alt.rd'),compress=T)

## 4. Check significance ##
apa.sig <- lapply(apa.alt,callsig)
#save(apa.sig,file=paste0(outdir,'apa.sig.rd'),compress=T)

legend_tab<-read.xlsx(paste0(inputdir,'legend.xlsx'))

for(i in 1:length(apa.sig))
{
tab<-apa.sig[[i]]$res
final_tab<-get_final_APA_tab(tab,comps$a[i],comps$b[i])

#sort
final_tab<-final_tab %>% arrange(.,peak_strand,peak_chr,peak_start)

legend<-gsub('group1',comps$a[i],legend_tab$legend)
legend<-gsub('group2',comps$b[i],legend)
legend<-cbind(colnames(final_tab),legend) %>% as.data.frame()
colnames(legend)<-c('columns','legend')

sig_tab<-filter_for_sig(final_tab)

tab_lst<-list(legend,final_tab,sig_tab)
names(tab_lst)<-c('legend','all_APA','sig_APA')
# Add legend
write.xlsx(tab_lst,paste0(outdir,'Differential_usage_w_sig_',compnames[i],'.xlsx'))


}



