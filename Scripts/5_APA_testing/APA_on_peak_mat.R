library(dplyr)
library(stringi)
library(stringr)
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
library(corrplot)

source('/home/sonas/Archive/scripts/final_scripts/scPASU_APA_testing_functions.R')


count_dir<-'/home/sonas/APA/output/4_APA/inputs/counts_mat_dir/'

f<-list.files(count_dir,full.names = TRUE)
peak_counts<-lapply(f,read.table,sep='\t', header=TRUE,stringsAsFactors=FALSE)


# Replace . with - in count mat 
for(i in 1:length(peak_counts))
{
  colnames(peak_counts[[i]])<-gsub('.','-',colnames(peak_counts[[i]]),fixed = TRUE)
  colnames(peak_counts[[i]])<-gsub('uro_','',colnames(peak_counts[[i]]))
}

# Merge counts
merged_counts<-do.call(cbind,peak_counts)

write.table(merged_counts,paste0('/home/sonas/APA/output/4_APA/inputs/u10_peak_counts.txt'),sep='\t')


# User inputs

#inputdir<-'/home/sonas/thesis_figures/ch4/APA/input/'
inputdir<-'/home/sonas/APA/output/4_APA/inputs/'

#outdir<-'/home/sonas/tingalab/Surbhi/PROJECTS_tinglab_drive/APA_Project/scPASU/output/4a_APA/outputs/'

#outdir<-'/home/sonas/thesis_figures/ch4/APA/test_output/'
outdir<-'/home/sonas/APA/output/4_APA/outputs/'

counts_file<-'u10_peak_counts.txt'

# Read merged counts
merged_counts<-read.table(paste0(inputdir,counts_file),sep='\t')
colnames(merged_counts)<-gsub('.','-',colnames(merged_counts),fixed=TRUE)

# Read per clus counts
# f<-list.files('/home/sonas/tingalab/Surbhi/PROJECTS_tinglab_drive/APA_Project/scPASU/output/3_PeakCounts/3c_CreatePeakMat_per_clus/counts_dir/',full.names = TRUE)
# 
# nm<-basename(f)
# nm<-gsub('2021_04_01_ureter10_uro_PC50_res0.2_','',nm)
# nm<-gsub('_peak_by_cell_count.rds','',nm)
# counts<-lapply(f,readRDS)


# Read meta data
meta<-read.xlsx(paste0(inputdir,'2021_04_01_ureter10_uro_PC50_res0.2_meta.xlsx'))


### Load results ###

## Peak ref is jtu$join and final_annotation is the peak column 
peak_ref_dir<-'/home/sonas/APA/output/2_PeakRef/tss/'
peak_file<-paste0(peak_ref_dir,'u10_uro_l300iA10mm3_e100_updated_tss_filt_final.txt')

peak_ref<-read.table(peak_file,header=TRUE)

# Replace peak name column with final annotation as this is more meaningful name 
col<-which(colnames(peak_ref)=='final_annotation')
peak_ref$peak<-peak_ref[,col]
peak_ref<-peak_ref[,-col]

# Remove P0 peaks

plist<-strsplit(peak_ref$peak,split=':')
tlist<-sapply(plist,'[',3)
rem<-which(tlist=='P0') # 4205 P0 peaks
peak_ref<-peak_ref[-rem,]
 

# Update peak ref to remove majority exon peak genes

# using bedops extract exon on HPC from https://bioinformatics.stackexchange.com/questions/18706/extract-only-exon-regions-from-gff-gtf-file-with-input-bed-regions

# gtf2bed < annotations.gtf | grep -wF exon > exons.bed

# OR use rtracklayer
library(rtracklayer)
library(Sierra)

gtf<-import('/home/sonas/APA/ref/genes.gtf')
exon   = subset(gtf, type=='exon')
peak_gr<-GRanges(peak_ref)

ov<-GenomicAlignments::findOverlaps(peak_gr, exon,type="within")
idx<-S4Vectors::queryHits(ov)

exon<-rep('other',nrow(peak_ref))
exon[idx]<-'exon'
peak_ref$exon<-exon


exon_count<-peak_ref %>% group_by(gene, exon) %>% tally()

# Retain only exon count
exon_count<-exon_count[which(exon_count$exon=='exon'),]

peak_count<-peak_ref %>% group_by(gene) %>% tally()
ecount<-rep(0,nrow(peak_count))
idx<-match(exon_count$gene,peak_count$gene)
ecount[idx]<-exon_count$n
peak_count$exon_count<-ecount
peak_count$exon_pct<-(peak_count$exon_count/peak_count$n)*100
peak_count<-peak_count[order(peak_count$exon_pct,decreasing = TRUE),]
  
write.xlsx(peak_count,'/home/sonas/APA/output/4_APA/outputs/ref_update/peak_exon_count.xlsx')

#### Archive ####
gtf_df<-as.data.frame(gtf)
exon<-gtf_df[which(gtf_df$type=='exon'),]

# Get chr,start,end,gene name, score and strand
exon_bed<-exon[,c(1,2,3,13,8,5)]
exon_bed$score<-rep(0,nrow(exon_bed))
write.table(x=exon_bed,file='exon.bed',col.names=FALSE,row.names=FALSE,quote=FALSE,sep='\t')

############# #######################


#rep_dir<-'/home/sonas/thesis_figures/ch4/data_partition/rep_dir/'

# Read sample and comparison tables #

comp_file<-paste0(inputdir,'comp_group.csv')
compdt <- fread(comp_file,header=TRUE)
comps <- data.table(a=compdt$group1,b=compdt$group2)
bc<-comps

# 5 at a time - run 26 and 28
comps<-comps[c(27,26),]

# Do it only for 2 comparison
#comps<-comps[c(14,15),]
# check for old comp 1st
#comps<-data.table("a"=c('b3','u4'),"b"=c('b2','b2'))
compnames <- paste0(comps$a,"_vs_",comps$b)

# Create stab

group<-c(comps$a,comps$b) %>% unique() %>% sort()
sample<-c(paste0(group,'_1'),paste0(group,'_2')) 
group<-rep(group,2)
n<-length(group)/2
stab<-cbind(sample,group) %>% as.data.frame()

#samp<-fread('/home/sonas/APA/scAPA/sample_sheet.csv')
#samp<-fread('/home/sonas/APA/scAPA/sample_sheet2.csv')

stab$group<- stab$group %>% as.factor()

cov<-10
#a<-'u4'
a<-comps$a
b<-comps$b
# a<-'b3'
# b<-'b2'
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
  ta <- testAPA2(stab=stab,meanmat=meanmat,peak_ref=peak_ref,peak_mat=peak_mat,a=a,b=b,adjust.var=NULL,ncpu=ncpu,min_peak=min_peak)
  
  # Run with batch correction
  # taadj <- testAPA2(stab=stab,meanmat=meanmat,peak_ref=peak_ref,peak_mat=peak_mat,a=a,b=b,adjust.var=adjust.var,ncpu=ncpu,min_peak=min_peak)
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



