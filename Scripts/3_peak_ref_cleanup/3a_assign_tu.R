#!/usr/bin/env Rscript

library(handy)
library(goldmine)
library(rtracklayer)
library(dplyr)
library(data.table)
library(gtools)

argv<-commandArgs(trailing = TRUE)
ncpu <- argv[1]
dir=argv[2] # path to scPASU dir where inputs and outputs are organized
peaks_dir=argv[3] # path to final peaks ref dir 
fprefix=argv[4] # file prefix
gtf_file=argv[5] # gtf file with full path

script_dir<-paste0(dir,'/scripts/')
source(paste0(script_dir,'/scPASU_functions.R'))

outdir<-paste0(dir,/output/3_RefinePeakRef/3a_assign_TU/')
inputdir<-paste0(dir,'/input/')

# HPC paths
#dir<-'/home/sonas/beegfs/APA/scPASU/input/'
#genedir<-'/home/sonas/beegfs/APA/scPASU/input/'
#genedir<-'/home/sonas/beegfs/APA/scAPA/ureter10/Archive/assign_tu'
#gtf_file<-"/home/sonas/beegfs/ref/refdata-gex-GRCh38-2020-A/genes/genes.gtf"


ensembl_biotype_table<-paste0(inputdir,"GRCh38_gene_biotypes.csv")
genes_file<-paste0(inputdir,"/genes.rds")
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
 genes<-readRDS(genes_file)
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

gtf<-gtf[!duplicated(gtf$gene_id),]
genes<-genes[genes$gene.id %in% gtf$gene_name,]

# Get gene type name from gtf
genes$biotype <- gtf[match(genes$gene.id,gtf$gene_name),]$gene_type


genes <- genes[biotype %in% bio.want,]
which((is.na(genes$biotype)) ==TRUE) %>% length()

red_ens <- reduceGenes(genes=genes,chrs=chrs,flank=5000,ncore=ncore)
save(red_ens,file=genes_file)
}

## Peaks file ##

peaks_dir_m=paste0(peaks_dir,'/minus/4_polya_supported_peak_ref/')
peaks_dir_p=paste0(peaks_dir,'/plus/4_polya_supported_peak_ref/')

files_m<-list.files(peaks_dir_m,full.names = TRUE)
file_lst_m<-lapply(files_m,fread,header=TRUE)
merged_m<-do.call(bind_rows,file_lst_m)

# correct colnames (if needed)
cols<-colnames(merged_m)

# In my outputs for minus strand end<start but this generates error in GRange so swap colnames
cat('Reordering columns in minus strand files \n')

#colnames(merged_m)<-cols[c(1,3,2,4:9,11,10,12:14)]

files_p<-list.files(peaks_dir_p,full.names = TRUE)
file_lst_p<-lapply(files_p,fread,header=TRUE)
merged_p<-do.call(bind_rows,file_lst_p)

merged<-rbind(merged_m,merged_p)

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

#merged<- read.table('/u10_uro_peak_universe.txt', sep='\t',header=TRUE)

cat('Creating GRange obj \n')
peaks<-makeGRanges(merged,strand=T)

cat('Assign TU \n')
jtu <- joinTus_peaks(allpeaks=peaks,rg=red_ens)

cat('Save all tables \n')
save(jtu,peaks,red_ens,file=output_file)


# Create Peak reference with all relevant columns merged after TU assignment #
cat('Add other relevant cols \n')

jtu$join<-as.data.frame(jtu$join)
# Remove multi-TU peaks
r1<-which(jtu$join$unique_peak==FALSE)
tu_tab<-jtu$join[-r1,] 
#r2<-which(jtu$join$unique_tu==FALSE)
#tu_tab<-jtu$join[-r2,] 

# Create peak per TU count
tu_peak<-jtu$join %>% dplyr::select(peak,tu) %>% group_by(tu) %>% tally()

# Peaks (saved from merging all peaks from MACS2 output)
peaks<-as.data.frame(jtu$polya_peaks)
found<-match(tu_tab$peak,peaks$peakID)
length(found)
matched_peaks<-peaks[found,]

# Sort all tables
tu_tab<-tu_tab[mixedorder(tu_tab$peak),]
matched_peaks<-matched_peaks[mixedorder(matched_peaks$peakID),]

# Now lets transfer peak info to TU table
identical(tu_tab$peak,matched_peaks$peakID) # TRUE

tu_tab$chr<-matched_peaks$seqnames
tu_tab$start<-matched_peaks$start
tu_tab$end<-matched_peaks$end
tu_tab$strand<-matched_peaks$strand
tu_tab$gene<-strsplit(tu_tab$tu_anno,split=':') %>% sapply(.,'[[',2)

# Also transfer other info for future use. This table can then serve as comprehensive features table

tu_tab$peak_width<-matched_peaks$peak_width
tu_tab$pr_chr<-matched_peaks$pr_chr
tu_tab$pr_start<-matched_peaks$pr_start
tu_tab$pr_end<-matched_peaks$pr_end
tu_tab$pr_strand<-matched_peaks$pr_strand
tu_tab$pr_width<-matched_peaks$pr_width
tu_tab$polya_count<-matched_peaks$polya_count
tu_tab$polya_count<-matched_peaks$peakID

# Now add peakID to tu annotation column (multi peak TU will look like TU1:gene:P1, TU1:gene:P2 and TU1:gene:P3 while single peak TU will look like TU2:gene:P0)

# Get TU count
tu_count<-tu_tab %>% group_by(tu) %>% tally()

# Sort
tu_count<-tu_count[order(tu_count$n),]

#Single peak per TU
single<-which(tu_count$n==1)
s<-match(tu_count$tu[single],tu_tab$tu)
single_tu<-tu_tab[s,]
single_tu$tu %>% unique() %>% length()
single_tu$final_annotation<-paste0(single_tu$tu_anno,':P0')

# Multi peak TU -
multi_tu<-tu_tab[-s,]

# Split by strand
multi_tu_p<-multi_tu[multi_tu$strand=='+',]
multi_tu_m<-multi_tu[multi_tu$strand=='-',]

# Now add final annotation col
multi_tu_p2<-create_final_annotation_col(multi_tu_p)
multi_tu_m2<-create_final_annotation_col(multi_tu_m,is_minus=TRUE)

# Merge plus and minus strand
multi_tu2<-rbind(multi_tu_p2,multi_tu_m2)

# Add P0 genes back to ref
merged_tu3<-rbind(single_tu,multi_tu2)

# Save
merged_tu3<-merged_tu3 %>% as.data.frame()

write.table(merged_tu3,paste0(outdir,fprefix,'_peak_universe_updated.txt'),sep='\t',row.names=FALSE,col.names=TRUE)




# Also save a strand specific .saf
saf_ref_m<-saf_ref[saf_ref$Strand=='-',]
saf_ref_p<-saf_ref[saf_ref$Strand=='+',]

write.table(saf_ref_m,paste0(outdir,fprefix,'_peak_universe_minus_updated.saf'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(saf_ref_p,paste0(outdir,fprefix,'_peak_universe_plus_updated.saf'),sep='\t',quote=FALSE,row.names=FALSE)



# Create SAF format peak ref too

# Select relevant columns
cols<-c('final_annotation','chr','start','end','strand')
select<-match(cols,colnames(merged_tu3))
saf_ref<-merged_tu3[,select]
colnames(saf_ref)<-c('GeneID','Chr','Start','End','Strand') 

cat('Creating SAF ref file \n')
write.table(saf_ref,paste0(outdir,fprefix,'_peak_universe_updated.saf'),sep='\t',quote=FALSE,row.names=FALSE)

