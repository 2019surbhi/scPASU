
#library(devtools)
#install_github("BMILAB/movAPA")

library(dplyr)
library(data.table)

library(reshape2)
library(RColorBrewer)
library(ggplot2)

library(ggbio)
library(IRanges)
library(GenomicRanges)
library(BSgenome)
library(rtracklayer)
library(GenomicFeatures)
library(DESeq2)
library(DEXSeq)
library(movAPA)


### Demo Data ###

#pac=read.csv('/Users/surbhisona/Downloads/arabidopsis_thaliana.SRP093950_amp.high_confidence.PAC.annotation.tpm.csv',stringsAsFactors =F)
# 
#pac <- pac[,-1]
# pac=dplyr::rename(pac, UPA_start = 'start', UPA_end='end', gene_type='biotype')
# 
# colData=as.data.frame(matrix(c('Amp','Amp','Amp'), ncol=1, dimnames=list(paste0('Amp311_R',1:3), 'group')))
# 
# PACds=readPACds(pacFile=pac, colDataFile=colData, noIntergenic=FALSE, PAname='PA')
# 
# # Try with subset - works!!
# pac_sub<-pac[,c(1:4,7,20:22)]
# PACds_sub=readPACds(pacFile=pac_sub, colDataFile=colData, noIntergenic=FALSE, PAname='PA')


## Next step ##

# library("BSgenome.Athaliana.TAIR.TAIR9")
# bsgenome <- BSgenome.Athaliana.TAIR.TAIR9
# 
# # Please make sure the chromosome name of your PAC data is the same as the BSgenome.
# seqnames(bsgenome) <- c(1:5,'Mt','Pt')
# seqnames(bsgenome)
# #> [1] "1"  "2"  "3"  "4"  "5"  "Mt" "Pt"
# 
# PACdsIP=removePACdsIP(PACds, bsgenome, returnBoth=TRUE,
#                       up=-10, dn=10, conA=6, sepA=7)
# 
# faFiles=faFromPACds(PACds, bsgenome, what='updn', fapre='updn',
#                     up=-300, dn=100, byGrp='ftr')
# faFiles=c("updn.3UTR.fa", "updn.CDS.fa", "updn.intergenic.fa", "updn.intron.fa")
# 
# # do it for entire ref
# faFile=faFromPACds(PACds, bsgenome, what='updn', fapre='updn',
#                     up=-300, dn=100)
# faFile<-"updn.fa"
# 
# ## Plot single nucleotide profiles using the extracted sequences and merge all plots into one.
# 
# plotATCGforFAfile ("updn.3UTR.fa", ofreq=FALSE, opdf=FALSE,
#                    refPos=301, mergePlots = TRUE)
# 
# 
# plotATCGforFAfile ("updn.fa", ofreq=FALSE, opdf=FALSE,
#                    refPos=301, mergePlots = TRUE)
# 

########## ms data ########################

# data("scPACds")
# head(scPACds@anno)


########### Preparing my data ################

library("BSgenome.Hsapiens.UCSC.hg38")
bsgenome <- BSgenome.Hsapiens.UCSC.hg38

# chr, UPA_start, UPA_end and strand - coordinates of peaks?
# coord - denotes anchored polyA site (set as peak end for + and peak start for -)
# few counts columns
ref_file<-'/Users/surbhisona/Desktop/movAPA/u10_uro_l300iA10mm3_e100_updated_tss_filt.txt'

ref<-read.table(ref_file,header=TRUE,stringsAsFactors = FALSE)

# subset relevant cols
ref2<-ref[,c("chr","start","end","strand","final_annotation","gene","peak_width","pr_start","pr_end","pr_width","polya_count")]

# rename

# Use peaks
#ref2<-dplyr::rename(ref2, UPA_start = 'start', UPA_end='end')


# Use PR - use this!
ref2<-dplyr::rename(ref2, 
                    UPA_start = 'pr_start', UPA_end='pr_end',
                    peak_start = 'start', peak_end='end')


# Add dummy counts for now
n<-nrow(ref2)
c1<-sample(1:100,n,replace = TRUE)
c2<-sample(1:100,n,replace = TRUE)
ref2$c1<-c1
ref2$c2<-c2

# Add coord

which(ref2$UPA_start>ref2$UPA_end)
ref2$coord<-ifelse(ref2$UPA_start<ref2$UPA_end,ref2$UPA_start,ref2$UPA_end)


# Prepare coldata

colData<-as.data.frame(matrix(c('dummy','dummy'), ncol=1, dimnames=list(paste0('c',1:2), 'group')))


# Create PACds
my_pac=readPACds(pacFile=ref2, colDataFile=colData, noIntergenic=FALSE, PAname='PA')


# Don't need this step : *removePACdsIP* since I already filtered internally primed reads

# Annotate

#gff=parseGff('/Users/surbhisona/Documents/PROJECTS/ref/refdata-gex-GRCh38-2020-A/genes/genes.gtf')

head(gff$anno.need)

my_pac=annotatePAC(my_pac, gff)

#faFiles=faFromPACds(my_pac, bsgenome, what='updn', fapre='updn',
#                    up=-200, dn=200)


### Stratified by annotation ####
faFiles=faFromPACds(my_pac, bsgenome, what='updn', fapre='updn',
                    up=-100, dn=100, byGrp='ftr')

faFiles=c("updn.3UTR.fa", "updn.CDS.fa", "updn.intergenic.fa", "updn.intron.fa")

# Plot single nucleotide profiles using the extracted sequences and merge all plots into one.

prefix<-'u10_uro_pr_'
pdf(paste0(prefix,'nucleotide_freq_by_group.pdf'),width=8,height=8)

plotATCGforFAfile (faFiles, ofreq=FALSE, opdf=FALSE,
                   refPos=101,mergePlots = TRUE)

dev.off()

#### All together #####
faFile=faFromPACds(my_pac, bsgenome, what='updn', fapre='updn',
                    up=-100, dn=100)
faFile<-"updn.fa"

#setwd('/Users/surbhisona/Desktop/movAPA/')
prefix<-'u10_uro_pr_'
pdf(paste0(prefix,'nucleotide_freq.pdf'),width=8,height=8)
plotATCGforFAfile (faFile, ofreq=FALSE, opdf=FALSE,
                   refPos=101,filepre = 'u10_uro')
dev.off()

## Plot single nucleotide profiles using the extracted sequences and merge all plots into one.
#plotATCGforFAfile (faFiles, ofreq=FALSE, opdf=FALSE,
#                   refPos=101, mergePlots = TRUE,)



#### polyA db #####

pdb_file<-'atlas.clusters.2.0.GRCh38.96.bed'
pdb<-read.table(pdb_file,sep='\t',header=FALSE)
colnames(pdb)<-c('chr','start','end','prID','avg_exp','strand','sample_prop','prot_count','avg_exp_all','anno','PAS')
head(pdb)

pdb$chr<-paste0('chr',pdb$chr)
chrs<-unique(gff$anno.frame$seqnames)
pdb<-pdb[pdb$chr %in% chrs,]

# Sanity check - end>start irrespective of strand
diff<-pdb$end-pdb$start
which(diff<0)

# Prepare data

pdb2<-dplyr::rename(pdb, 
                    UPA_start = 'start', UPA_end='end')


# Add dummy counts for now
n<-nrow(pdb2)
c1<-sample(1:100,n,replace = TRUE)
c2<-sample(1:100,n,replace = TRUE)
pdb2$c1<-c1
pdb2$c2<-c2

# Add coord

# coord col 
pa<-strsplit(pdb2$prID,split=':') %>% sapply(.,'[[',2) %>% as.numeric()
pa<-pa-1
#pdb2$coord<-pa

# Use start (similar to my data)
pdb2$coord<-ifelse(pdb2$UPA_start<pdb2$UPA_end,pdb2$UPA_start,pdb2$UPA_end)


# Prepare coldata

colData<-as.data.frame(matrix(c('dummy','dummy'), ncol=1, dimnames=list(paste0('c',1:2), 'group')))


# Create PACds
pdb_pac=readPACds(pacFile=pdb2, colDataFile=colData, noIntergenic=FALSE, PAname='PA')

pdb_pac=annotatePAC(pdb_pac, gff)

#faFiles=faFromPACds(my_pac, bsgenome, what='updn', fapre='updn',
#                    up=-200, dn=200)


### Stratified by annotation ####
faFiles=faFromPACds(pdb_pac, bsgenome, what='updn', fapre='updn',
                    up=-100, dn=100, byGrp='ftr')

faFiles=c("updn.3UTR.fa", "updn.CDS.fa", "updn.intergenic.fa", "updn.intron.fa")

# Plot single nucleotide profiles using the extracted sequences and merge all plots into one.

prefix<-'pdb_hs_hg38_'
prefix<-'pdb_hs_hg38_using_start_'
pdf(paste0(prefix,'nucleotide_freq_by_group.pdf'),width=8,height=8)
plotATCGforFAfile (faFiles, ofreq=FALSE, opdf=FALSE,
                   refPos=101,mergePlots = TRUE)

dev.off()

#### All together #####
faFile=faFromPACds(pdb_pac, bsgenome, what='updn', fapre='updn',
                   up=-100, dn=100)
faFile<-"updn.fa"

setwd('/Users/surbhisona/Desktop/movAPA/')
prefix<-'pdb_hs_hg38_'
prefix<-'pdb_hs_hg38_using_start_'

pdf(paste0(prefix,'nucleotide_freq.pdf'),width=8,height=8)
plotATCGforFAfile (faFile, ofreq=FALSE, opdf=FALSE,
                   refPos=101,filepre = 'u10_uro')
dev.off()




tab<-fread('/Users/surbhisona/Downloads/scPASU ref/2021_01_18_ureter10_sierra_merged_peak_annotations.txt')
colnames(tab)[1]<-'pA'
head(tab)


sierra<-dplyr::rename(tab, 
                    UPA_start = 'start', UPA_end='end',chr='seqnames')


# Add dummy counts for now
n<-nrow(sierra)
c1<-sample(1:100,n,replace = TRUE)
c2<-sample(1:100,n,replace = TRUE)
sierra$c1<-c1
sierra$c2<-c2

# Add coord

which(sierra$UPA_start>sierra$UPA_end)
sierra$coord<-ifelse(sierra$UPA_start<sierra$UPA_end,sierra$UPA_start,sierra$UPA_end)


# Prepare coldata


colData_sierra<-as.data.frame(matrix(c('dummy','dummy'), ncol=1, dimnames=list(paste0('c',1:2), 'group')))


sierra<-as.data.frame(sierra)

sierra_pac<-readPACds(pacFile=sierra, colDataFile=colData_sierra, noIntergenic=FALSE, PAname='PA')

faFile_sierra=faFromPACds(sierra_pac, bsgenome, what='updn', fapre='updn',
                   up=-100, dn=100)
faFile<-"updn.fa"

setwd('/Users/surbhisona/Desktop/movAPA/')
prefix<-'pdb_sierra_'
prefix<-'pdb_sierra_using_start_'

pdf(paste0(prefix,'nucleotide_freq.pdf'),width=8,height=8)
plotATCGforFAfile (faFile, ofreq=FALSE, opdf=FALSE,
                   refPos=101,filepre = 'u10_uro_sierra')
dev.off()



################### PR by annotation ##################

pdb_pac
my_pac

# Annotate Sierra PAC as well
sierra_pac=annotatePAC(sierra_pac, gff)

# Get stats
pstats_pdb=movStat(pdb_pac, minPAT=c(1, 5, 10, 20, 50, 60), ofilePrefix=NULL)
pstats_mypac=movStat(my_pac, minPAT=c(1, 5, 10, 20, 50, 60), ofilePrefix=NULL)
pstats_sierra=movStat(sierra_pac, minPAT=c(1, 5, 10, 20, 50, 60), ofilePrefix=NULL)

plotPACdsStat(pstats_pdb, pdfFile='PACds_stat_total_pdb.pdf', 
             minPAT=c(5,10), conds=c('total'))

plotPACdsStat(pstats_mypac, pdfFile='PACds_stat_total_scPASU.pdf',
              minPAT=c(5,10),conds=c('total'))

plotPACdsStat(pstats_sierra, pdfFile='PACds_stat_total_sierra.pdf',
              minPAT=c(5,10),conds=c('total'))

