#!/usr/bin/env Rscript

library(dplyr)
library(data.table)

args <- commandArgs(trailing = TRUE)
dir=args[1]
strand=args[2]
chr=args[3]

#dir<-'/home/sonas/beegfs/results/scAPA/my_scAPA/ureter_uro/macs2_u10_uro_merged/pr_table_minus/'
#dir<-'/home/sonas/beegfs/results/scAPA/my_scAPA/ureter_uro/macs2_u10_uro_merged/pr_table_minus_genomicA_filtered/'
#strand="minus"

#dir<-'/home/sonas/beegfs/results/scAPA/my_scAPA/ureter_uro/macs2_u10_uro_merged/pr_table_plus/'
#dir<-'/home/sonas/beegfs/results/scAPA/my_scAPA/ureter_uro/macs2_u10_uro_merged/pr_table_plus_genomicA_filtered/'
#strand="plus"

out<-paste0(dir,'final_ref_table/')
prdir<-paste0(dir,'pr_dir/')
pr_tab_dir<-paste0(dir,'/pr_table_dir/')
count_dir<-paste0(dir,'peaks_count_dir/')

pr_file<-paste0(prdir,'pr_',chr,'_sorted_',strand,'.bed')
pr_tab_file<-paste0(pr_tab_dir,'pr_table_',chr,'_',strand,'.bed')
count_file<-paste0(count_dir,'peaks_',chr,'_count_sorted_',strand,'.txt')

# Occassionally PR tables are empty
#if(file.size(pr_file)==0)
#{
#  cat('Empty file ')
#  quit(save = 'no')
#}



pr<-read.table(pr_file,header=FALSE,sep="\t")
if(strand=='minus')
{ colnames(pr)<-c('peakID','end','start')
 }else if(strand=='plus')
{
  colnames(pr)<-c('peakID','start','end')
}

pr_tab<-read.table(pr_tab_file,header=FALSE,sep="\t")

if(strand=='minus')
  {colnames(pr_tab)<-c('peak_chr','peak_end','peak_start','peakID','peak_strand','pr_chr','polya_end','polya_start','readID','A_len','pr_strand','polya')
}else if(strand=='plus'){
  colnames(pr_tab)<-c('peak_chr','peak_start','peak_end','peakID','peak_strand','pr_chr','polya_start','polya_end','readID','A_len','pr_strand','polya')
}else{
  cat('Error: must specify strand \n')
}
 
count<-fread(count_file)
colnames(count)<-c('count','peakID')
 
# Remove polya coordinate columns
pr_tab<-pr_tab[,-c(7,8,9,10)]
 
pr_tab$pr_end<-0
pr_tab$pr_start<-0
 
indx<-match(pr$peakID,pr_tab$peakID)
pr_tab$pr_end[indx]<-pr$end
pr_tab$pr_start[indx]<-pr$start

if(strand=='minus')
{
pr_tab$peak_width=pr_tab$peak_start-pr_tab$peak_end
pr_tab$pr_width=pr_tab$pr_start-pr_tab$pr_end
}else if(strand=='plus')
{
pr_tab$peak_width=pr_tab$peak_end-pr_tab$peak_start
pr_tab$pr_width=pr_tab$pr_end-pr_tab$pr_start
}

#Add polya count
pr_tab$polya_count<-0
indx2<-match(count$peakID,pr_tab$peakID)
pr_tab$polya_count[indx2]<-count$count
 
#append final table
#final_tab<-rbind(final_tab,pr_tab)
 
#save individual chr output
fname<-basename(pr_tab_file) %>% gsub('pr_table','ref_table_updated',.)
 
write.table(pr_tab,file = paste0(out,fname),row.names = FALSE,sep = "\t")

