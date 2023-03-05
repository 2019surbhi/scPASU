#!/usr/bin/env Rscript

library(dplyr)
library(data.table)
library(gtools)

argv<-commandArgs(trailing = TRUE)
peak_ref_dir=argv[1] # path to final peak ref table (after assinging TU)
outdir=argv[2] # output dir for final peak ref
fprefix=argv[3] # file prefix


# Example paths on HPC
#outdir<-'/home/sonas/beegfs/APA/scAPA/ureter10/peak_mat/final_ref/'


#peak_ref_dir<-'/home/sonas/beegfs/APA/scAPA/ureter10/assign_tu/u10_uro_unfiltered_3polya/'
#fprefix<-'u10_uro_unfiltered_3polya'


# Create dir if doesn't exist

if(!dir.exists(outdir))
        {dir.create(outdir)}

# Load data
load(paste0(peak_ref_dir,fprefix,'_jtu.rds'))

ngenes<-jtu$join$tu %>% unique() %>% length() # 13114 unique TU
npeaks<-jtu$join$peak %>% unique() %>% length() # 89048 unique peaks

# Single peak TUs
tu_peak<-jtu$join %>% select(peak,tu) %>% group_by(tu) %>% tally()
single_peak_tu<-which(tu_peak$n==1) %>% length() # 3295

# Multi peak TU
multi_peak_tu<-which(tu_peak$n>1) %>% length() # 9819


# Further breakdown this count
multi_peak_tu_10<-which(tu_peak$n>1 & tu_peak$n<=10) %>% length() # 7149
multi_peak_tu_20<-which(tu_peak$n>10 & tu_peak$n<=20) %>% length()  # 447
multi_peak_tu20plus<-which(tu_peak$n>20) %>% length() # 136


# Remove multi-TU peaks
rem<-which(jtu$join$unique_peak==FALSE)
tu_tab<-jtu$join[-rem,]

# TU after filtering
ngenes2<-tu_tab$tu %>% unique() %>% length() # 12379 unique TU
npeaks2<-tu_tab$peak %>% unique() %>% length() # 87226 unique peaks

# Single peak TUs after filtering
tu_peak2<-tu_tab %>% select(peak,tu) %>% group_by(tu) %>% tally()
single_peak_tu2<-which(tu_peak2$n==1) %>% length() # 3195

# Multi peak TU
multi_peak_tu2<-which(tu_peak2$n>1) %>% length() # 9184


# Further breakdown this count
multi_peak_tu_10f<-which(tu_peak2$n>1 & tu_peak2$n<=10) %>% length() # 7149
multi_peak_tu_20f<-which(tu_peak2$n>10 & tu_peak2$n<=20) %>% length()  # 447
multi_peak_tu20plusf<-which(tu_peak2$n>20) %>% length() # 136

stats<-data.frame('Peaks'=npeaks,'TU'=ngenes,'SinglePeakTU'=single_peak_tu,'MultiPeakTU'=multi_peak_tu,
'PeaksTU_2to10'=multi_peak_tu_10,'PeaksTU_11to20'=multi_peak_tu_20,'PeaksTU>20'=multi_peak_tu20plus)

filtered_stats<-data.frame('Peaks'=npeaks2,'TU'=ngenes2,'SinglePeakTU'=single_peak_tu2,'MultiPeakTU'=multi_peak_tu2,
'PeaksTU_2to10'=multi_peak_tu_10f,'PeaksTU_11to20'=multi_peak_tu_20f,'PeaksTU>20'=multi_peak_tu20plusf)

final_stats<-rbind(stats,filtered_stats)
 
rownames(final_stats)<-c('with_all_peaks','ambiguous_peaks_removed')
# save
write.table(final_stats,paste0(outdir,fprefix,'_stats.txt'),sep='\t',row.names=TRUE,col.names=TRUE)


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


# Also transfer other info for future use. This table can then serve as comprehensive features table

tu_tab$peak_width<-matched_peaks$peak_width
tu_tab$pr_chr<-matched_peaks$pr_chr
tu_tab$pr_start<-matched_peaks$pr_start
tu_tab$pr_end<-matched_peaks$pr_end
tu_tab$pr_strand<-matched_peaks$pr_strand
tu_tab$pr_width<-matched_peaks$pr_width
tu_tab$polya_count<-matched_peaks$polya_count

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
single_tu$tu_anno2<-paste0(single_tu$tu_anno,':P0')

# Multi peak TU -
multi_tu<-tu_tab[-s,]

# Split by strand
multi_tu_p<-multi_tu[multi_tu$strand=='+',]
multi_tu_m<-multi_tu[multi_tu$strand=='-',]

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

multi_tu_p<-create_final_annotation_col(multi_tu_p)
multi_tu_m<-create_final_annotation_col(multi_tu_m)

multi_tu<-rbind(multi_tu_p,multi_tu_m)

colnames(single_tu)<-gsub('tu_anno2','final_annotation',colnames(single_tu))
merged_tu<-rbind(single_tu,multi_tu)

# Save
merged_tu<-merged_tu %>% as.data.frame()

write.table(merged_tu,paste0(outdir,fprefix,'_updated.txt'),sep='\t',row.names=FALSE,col.names=TRUE)
