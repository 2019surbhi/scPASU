library(rtracklayer)
library(dplyr)
library(GenomicFeatures)
library(GenomicDistributions)
library(data.table)
library(ggplot2)
library(readxl)
library(openxlsx)


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


create_final_annotation_col<-function(multi_tu)
{
  # sort
  multi_tu<-multi_tu[mixedorder(multi_tu$tu),]
  
  # Sort minus strand differently for correct annotation
  if(multi_tu$strand=='-')
  {
    multi_tu<-multi_tu %>% arrange(., chr,desc(start))
  }
  
  tu<-unique(multi_tu$tu)
  cnt<-multi_tu %>% group_by(tu) %>% tally()
  cnt<-cnt[order(cnt$n),]
  
  multi_tu$final_annotation<-rep('none',nrow(multi_tu))
  
  for(i in 1:length(tu))
  {
    
    indx<-match(tu[i],cnt$tu)
    
    n<-cnt$n[indx]
    p<-paste0('P',(1:n))
    indx2<-which(multi_tu$tu %in% tu[i])
    new_anno<-paste0(multi_tu$tu_anno[indx2],':',p)
    multi_tu$final_annotation[indx2]<-new_anno
    
  }
  return(multi_tu)
  
}



load('/home/sonas/tingalab/Surbhi/PROJECTS_tinglab_drive/APA_Project/scPASU/output/2_PeakRef/2d_TU_annotated/genes.rds')

extn<-100

out<-'/home/sonas/APA/tss/'
fprefix<-paste0('u10_uro_l300iA10mm3_e',extn)

peak_ref<-read.table(paste0('/home/sonas/tingalab/Surbhi/PROJECTS_tinglab_drive/APA_Project/scPASU/output/2_PeakRef/2d_TU_annotated/',fprefix,'_updated.txt'),header=TRUE)

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

# Create gene size
new_ref$gene_size<-ifelse(new_ref$strand=='-',(new_ref$tss-new_ref$tu_end), (new_ref$tu_end-new_ref$tss))


# Calculate TSS distance
new_ref$tss_dist<-ifelse(new_ref$strand=='-',(new_ref$tss-new_ref$start), (new_ref$start-new_ref$tss))

# Calculate 20% gene size 
new_ref$gene_size_20p<-round((new_ref$gene_size*0.2),0)


# Now filter
for(i in 1:nrow(new_ref))
{
  if(new_ref$gene_size_20p[i]<=1000)
  {
    if(new_ref$tss_dist[i]<=new_ref$gene_size_20p[i])
    {
      new_ref$filter_peak[i]<-TRUE
    }else{
      new_ref$filter_peak[i]<-FALSE
    }
  }else if(new_ref$tss_dist[i]<=1000)
  {
    new_ref$filter_peak[i]<-TRUE
  }else{
    new_ref$filter_peak[i]<-FALSE
  }
}

filt_ref<-new_ref[new_ref$filter_peak==FALSE,]

# Update peaks

peak_count<-filt_ref %>% group_by(tu) %>% tally()

# Sort
peak_count<-peak_count[order(peak_count$n),]

#Single peak per TU
single<-which(peak_count$n==1)
s<-match(peak_count$tu[single],filt_ref$tu)
filt_ref_p0<-filt_ref[s,]

# Extract multi peak TU
filt_ref_pn<-filt_ref[-s,]

# Split by strand
filt_ref_pn_p<-filt_ref_pn[filt_ref_pn$strand=='+',]
filt_ref_pn_m<-filt_ref_pn[filt_ref_pn$strand=='-',]

filt_ref_pn_p2<-create_final_annotation_col(filt_ref_pn_p)
filt_ref_pn_m2<-create_final_annotation_col(filt_ref_pn_m)

filt_ref_pn2<-rbind(filt_ref_pn_p2,filt_ref_pn_m2)

final_ref<-rbind(filt_ref_p0,filt_ref_pn2)


# Correct minus strand annotation

start<-ifelse(final_ref$strand=='+',final_ref$start,final_ref$end)
end<-ifelse(final_ref$strand=='+',final_ref$end,final_ref$start)

final_ref$start<-start
final_ref$end<-end


write.table(final_ref,paste0(out,fprefix,'_updated_tss_filt_final.txt'),sep='\t',row.names=FALSE)







