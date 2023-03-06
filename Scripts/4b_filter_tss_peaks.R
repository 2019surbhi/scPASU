library(dplyr)
library(openxlsx)


## function ##


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


create_peak_ref_annotation<-function(tu_tab)
{
tu_count<-tu_tab %>% group_by(tu) %>% tally()

# Sort
tu_count<-tu_count[order(tu_count$n),]

#Single peak per TU
single<-which(tu_count$n==1)
s<-match(tu_count$tu[single],tu_tab$tu)
single_tu<-tu_tab[s,]
p0_count<-single_tu$tu %>% unique() %>% length()
single_tu$tu_anno2<-paste0(single_tu$tu_anno,':P0')

single_tu$final_annotation<-single_tu$tu_anno2
rem_col<-which(colnames(single_tu)=='tu_anno2')
single_tu<-single_tu[,-rem_col]


# Multi peak TU -
multi_tu<-tu_tab[-s,]

# Split by strand
multi_tu_p<-multi_tu[multi_tu$strand=='+',]
multi_tu_m<-multi_tu[multi_tu$strand=='-',]


multi_tu_p<-create_final_annotation_col(multi_tu_p) 
multi_tu_m<-create_final_annotation_col(multi_tu_m)

multi_tu<-rbind(multi_tu_p,multi_tu_m)

merged_tu<-rbind(single_tu,multi_tu)

merged_tu<-merged_tu %>% as.data.frame()

if('peak' %in% colnames(merged_tu))
{
  merged_tu$peak<-strsplit(merged_tu$final_annotation,split=':') %>%
    sapply(.,'[[',3)
}

return(list(merged_tu,p0_count))

}


# Read peak ref

#peak_file<-'/home/sonas/tingalab/Surbhi/PROJECTS_tinglab_drive/APA_Project/Ureter_uro/plot and stats/u10_uro_unfiltered_3polya_updated_new.txt'

extn<-25
extn<-50
extn<-75
extn<-100
extn<-125
extn<-149

out<-'/home/sonas/APA/tss/'
fprefix<-paste0('u10_uro_l300iA10mm3_e',extn)

peak_file<-paste0('/home/sonas/tingalab/Surbhi/PROJECTS_tinglab_drive/APA_Project/scPASU/output/2_PeakRef/2d_TU_annotated/',fprefix,'_updated.txt')

peak_ref<-read.table(peak_file,header=TRUE)

peak_ref$peak<-strsplit(peak_ref$final_annotation,split=':') %>% sapply(.,'[[',3)
peak_ref$gene<-strsplit(peak_ref$final_annotation,split=':') %>% sapply(.,'[[',2)


# Read TSS table

tss_tab<-read.xlsx(paste0(out,fprefix,'_tss_dist.xlsx'))
idx<-which(tss_tab$filter_p1==1)
filt_tss<-tss_tab[idx,c('chr','start','end','strand','gene')]

peak_ref<-peak_ref %>% arrange(.,strand,chr,start)
filt_tss<-filt_tss %>% arrange(.,strand,chr,start)

# Find genes whose P1 peaks need to be filtered
int<-intersect(filt_tss$gene,peak_ref$gene)
pr_keep<-peak_ref[!(peak_ref$gene %in% int),]
pr_rem<-peak_ref[(peak_ref$gene %in% int),]

rem<-which(pr_rem$peak=='P1')
pr_rem<-pr_rem[-rem,]

# Now update annotation
pr_rem %>% head()
lst<-create_peak_ref_annotation(pr_rem) 
pr_rem2<-lst[[1]]
p0_cnt<-lst[[2]]
peak_ref2<-rbind(pr_keep,pr_rem2)


write.table(peak_ref2,paste0(out,fprefix,'_updated_tss_filt.txt'),sep='\t',row.names=FALSE)

