library(rtracklayer)
library(dplyr)
library(GenomicFeatures)
library(GenomicDistributions)
library(data.table)
library(ggplot2)
library(readxl)
library(openxlsx)



##### Get TSS from TU ####

# Swap strand
fix_start<-function(tab)
{
  minus_idx<-which(tab$strand=='-')
  start<-tab$start[minus_idx]
  end<-tab$end[minus_idx]
  
  tab$start[minus_idx]<-end
  tab$end[minus_idx]<-start
  
  
  return(tab)
}




load('/home/sonas/tingalab/Surbhi/PROJECTS_tinglab_drive/APA_Project/scPASU/output/2_PeakRef/2d_TU_annotated/genes.rds')

peak_ref<-read.table('/home/sonas/tingalab/Surbhi/PROJECTS_tinglab_drive/APA_Project/scPASU/output/2_PeakRef/2d_TU_annotated/u10_uro_l300iA10mm3_e100_updated.txt',header=TRUE)

out<-'/home/sonas/APA/tss/'
fprefix<-'u10_uro_l300iA10mm3_e100'



new_ref  <- data.table(peak_ref, key="tu")

# Add relevant columns by TU #
new_ref <- new_ref[,list(peak=peak,
                         coding=coding,
                         unique_peak=unique_peak,
                         unique_tu=unique_tu,
                         num_peak=length(peak),
                         chr=chr, start=start, end=end,
                         strand=strand, final_annotation=final_annotation),
                   by="tu"]

# Remove redundant rows
genes<-strsplit(new_ref$final_annotation,split=':') %>% sapply(.,'[[',2) 
new_ref$gene<-genes
uniq<-unique(new_ref$gene)
new_ref<-new_ref[match(uniq,new_ref$gene),]

# Fix start for minus strand
new_ref<-fix_start(new_ref)

# check if it worked
which(new_ref$start[new_ref$strand=='-']< new_ref$end[new_ref$strand=='-']) %>% length() 

# Read TU table
tu<-red_ens$tu %>% as.data.frame()
tu<-fix_start(tu)

# check if it worked
which(tu$start[tu$strand=='-']< tu$end[tu$strand=='-']) %>% length() 

# Add TSS from TU table
idx<-match(new_ref$tu,tu$tu)
new_ref$tss<-tu$start[idx]
new_ref$tu_end<-tu$end[idx]

new_ref_p0 <- new_ref[num_peak==1,c('chr','start','end','strand','final_annotation','gene','tss','tu_end')] 
new_ref_pn <- new_ref[num_peak>1,c('chr','start','end','strand','final_annotation','gene','tss','tu_end')] 

new_ref_pn$peak<-strsplit(new_ref_pn$final_annotation,split=':') %>% sapply(.,'[[',3)

# Save only P1 peaks
new_ref_pn_p1<-new_ref_pn[new_ref_pn$peak=='P1',]
#tss_cleaned_p<-tss_cleaned_p[tss_cleaned_p$gene %in% new_ref_pn_p1$gene,]

# Create TSS table
tss_p1<-new_ref_pn_p1[,c('chr','tss','tu_end','strand','gene','start','end')]


# Split by strand 
tss_p1_p<-tss_p1[tss_p1$strand=='+']
tss_p1_m<-tss_p1[tss_p1$strand=='-']

# Create gene size
tss_p1_p$gene_size<-tss_p1_p$tu_end-tss_p1_p$tss
tss_p1_m$gene_size<-tss_p1_m$tss-tss_p1_m$tu_end

# Calculate tss distance
tss_dist_p<-tss_p1_p$start-tss_p1_p$tss
tss_p1_p$tss_dist<-tss_dist_p

tss_dist_m<-tss_p1_m$tss-tss_p1_m$start
tss_p1_m$tss_dist<-tss_dist_m

tss_p1<-rbind(tss_p1_p,tss_p1_m)

# 5% of gene size
tss_p1$gene_5p<-(tss_p1$gene_size * 0.05)

# Distance less than 500 from tss - mark for filtering
idx<-which(tss_p1$tss_dist<=500)
tss_p1$tss_dist_less_than_500bp<-rep(0,nrow(tss_p1))
tss_p1$tss_dist_less_than_500bp[idx]<-1

# 5p gene size <=500 
idx2<-which(tss_p1$gene_5p<=500)
tss_p1$gene_size5p_less_than500<-rep(0,nrow(tss_p1))
tss_p1$gene_size5p_less_than500[idx2]<-1

# 5p gene size <=500 and tss dist <5p gene size
idx3<-which((tss_p1$gene_5p<=500) & (tss_p1$tss_dist<=tss_p1$gene_5p))
tss_p1$tss_dist_less_than_5p_gene_size<-rep(0,nrow(tss_p1))
tss_p1$tss_dist_less_than_5p_gene_size[idx3]<-1

filter_p1<-rep(0,nrow(tss_p1))
# Now filter
for(i in 1:nrow(tss_p1))
{
  if(tss_p1$gene_5p[i]<=500)
  {
    if(tss_p1$tss_dist[i]<=tss_p1$gene_5p[i])
    {
      filter_p1[i]<-1
    }
  }else if(tss_p1$tss_dist[i]<=500)
    {
    filter_p1[i]<-1
  }
}

tss_p1$filter_p1<-filter_p1

tss_p1<-tss_p1 %>% arrange(chr,tss)

# Sanity check
# tss_p1<-tss_p1 %>% arrange(gene)
# new_ref_pn_p1<-new_ref_pn_p1 %>% arrange(gene)
# identical(tss_p1$gene,new_ref_pn_p1$gene)
# 
# identical(tss_p1$gene[tss_p1$strand=='+'],new_ref_pn_p1$gene[new_ref_pn_p1$strand=='+'])
# 
# identical(tss_p1$gene[tss_p1$strand=='-'],new_ref_pn_p1$gene[new_ref_pn_p1$strand=='-'])



tss_p1<-tss_p1 %>% arrange(chr,tss_dist)
write.xlsx(x = tss_p1,paste0(out,fprefix,'tss_dist.xlsx'))


