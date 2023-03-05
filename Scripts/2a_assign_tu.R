#!/usr/bin/env Rscript

library(handy)
library(goldmine)
library(rtracklayer)
library(dplyr)
library(data.table)

#source('/home/sonas/beegfs/APA/scAPA/ureter10/misc/scAPA_paddle.R')
#source('/home/sonas/beegfs/results/scAPA/my_scAPA/ureter_uro/macs2_u10_uro_merged/script/scAPA_paddle.R')

argv<-commandArgs(trailing = TRUE)
ncpu <- argv[1]
dir=argv[2] # path to scPASU dir where inputs and outputs are organized
peaks_dir=argv[3] # path to final peaks ref dir 
fprefix=argv[4] # file prefix
gtf_file=argv[5] # gtf file with full path

script_dir<-paste0(dir,'/scripts/')
source(paste0(script_dir,'/scPASU_functions.R'))

outdir<-paste0(dir,/output/2_PeakRef/2d_TU_annotated/')
inputdir<-paste0(dir,'/input/')
# HPC paths
#dir<-'/home/sonas/beegfs/APA/scPASU/input/'
#genedir<-'/home/sonas/beegfs/APA/scPASU/input/'
#genedir<-'/home/sonas/beegfs/APA/scAPA/ureter10/Archive/assign_tu'
#gtf_file<-"/home/sonas/beegfs/ref/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
ensembl_biotype_table<-paste0(inputdir,"GRCh38_gene_biotypes.csv")
genes_file<-paste0(outdir,"/genes.rds")
output_file<-paste0(outdir,"/",fprefix,"_jtu.rds")
track_dir<-paste0(outdir,"/tracks/")
cache_dir<-paste0(outdir,"/tmp/")
#input_file<-paste0(dir,"/bam_path.csv")
chrs <- handy::chrs()
dist <- 10
cov <- 10
# Create dir if doesn't exist

if(!dir.exists(outdir))
        {dir.create(outdir)}

if(!dir.exists(track_dir))
        {dir.create(track_dir)}

if(!dir.exists(cache_dir))
        {dir.create(cache_dir)}

### Create genes table ###
if(file.exists(genes_file))
{
 load(genes_file)
}else{
genes <- getGenes2("gencode",genome="hg38", gencodetable="wgEncodeGencodeCompV32", cachedir=cache_dir,sync=FALSE,url='https://hgdownload.soe.ucsc.edu/goldenPath/')

#dim(genes)
#head(genes)

# Exclude genes in scaffold
genes <- genes[chr %in% chrs,]

#dim(genes)

# Read biotable
tab <- fread(ensembl_biotype_table)

# Select gene types for coding and long non-coding
bio.want <- tab[tab$group %in% c("Protein_coding","Long_noncoding"),]$Var1

g<-rtracklayer::import(gtf_file)
gtf=as.data.frame(g)

gtf<-gtf[!duplicated(gtf$gene_name),]
genes<-genes[genes$gene.id %in% gtf$gene_name,]

# Get gene type name from gtf
genes$biotype <- gtf[match(genes$gene.id,gtf$gene_name),]$gene_type


genes <- genes[biotype %in% bio.want,]
which((is.na(genes$biotype)) ==TRUE) %>% length()

red_ens <- reduceGenes(genes=genes,chrs=chrs,flank=5000,ncore=ncore)
save(red_ens,file=genes_file)
}

## Peaks file ##

peaks_dir_m=paste0(peaks_dir,'minus/filtered_ref_table/')
peaks_dir_p=paste0(peaks_dir,'plus/filtered_ref_table/')

files_m<-list.files(peaks_dir_m,full.names = TRUE)
selected<-files_m[grep('chr',files_m)]
file_lst_m<-lapply(selected,read.table,sep='\t',header=TRUE)
merged_m<-do.call(bind_rows,file_lst_m)

# correct colnames (if needed)
cols<-colnames(merged_m)

# In my outputs for minus strand end<start but this generates error in GRange so swap colnames
cat('Reordering columns in minus strand files \n')
colnames(merged_m)<-cols[c(1,3,2,4:8,10,9,11:13)]

files_p<-list.files(peaks_dir_p,full.names = TRUE)
selected<-files_p[grep('chr',files_m)]
file_lst_p<-lapply(selected,read.table,sep='\t',header=TRUE)
merged_p<-do.call(bind_rows,file_lst_p)

merged<-bind_rows(merged_m,merged_p)

# Rename columns since GRange recognizes start, end and strand only
cat('Renaming colnames to create GRange obj \n')

colnames(merged)<-gsub('\\bpeak_chr\\b','chr',colnames(merged))
colnames(merged)<-gsub('\\bpeak_start\\b','start',colnames(merged))
colnames(merged)<-gsub('\\bpeak_end\\b','end',colnames(merged))
colnames(merged)<-gsub('\\bpeak_strand\\b','strand',colnames(merged))

# Also fix peak name i.e. remove *_plus_ prefix (not needed if MACS output are correctly generated)
#merged$peakID<-gsub('\\b*_plus_','',merged$peakID)
#merged$peakID<-gsub('[*]','',merged$peakID)

# Save the dataframe
cat('Saving table \n')
write.table(merged,paste0(outdir,fprefix,'_peak_universe.txt'),sep='\t',row.names=FALSE,col.names=TRUE)

#convert to GRange

# Read in the data file

#merged<- read.table('/home/sonas/beegfs/APA/scAPA/ureter10/u10_uro_BAMfiltered_peak_universe/u10_uro_BAMfiltered_peak_universe.txt', sep='\t',header=TRUE)

cat('Creating GRange obj \n')
peaks<-makeGRanges(merged,strand=T)

cat('Assign TU \n')
jtu <- joinTus_peaks(allpeaks=peaks,rg=red_ens)

cat('Save all tables \n')
save(jtu,peaks,red_ens,file=output_file)
