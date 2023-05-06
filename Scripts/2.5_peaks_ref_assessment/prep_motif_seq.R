library(dplyr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Rsamtools)
library(Biostrings)
library(ShortRead)
#library(data.table)

ref_file<-'/Volumes/tingalab/Surbhi/PROJECTS_tinglab_drive/APA_Project/scPASU/output/2_PeakRef/2e_TSS_filtered/u10_uro_l300iA10mm3_e100_updated_tss_filt.txt'

ref<-read.table(ref_file,header=TRUE,stringsAsFactors = FALSE)

amb<-which(ref$strand!=ref$pr_strand)
ref$pr_strand[amb]<-ref$strand[amb]

pr_ref<-dplyr::rename(ref,
                      peak_start='start',
                      peak_end='end',
                      peak_strand='strand')

pr_ref<-dplyr::rename(pr_ref,
                      start='pr_start',
                      end='pr_end',
                      strand='pr_strand')
pr_gr<-makeGRangesFromDataFrame(pr_ref,keep.extra.columns = TRUE)

us<-50
pr_gr2<-resize(pr_gr,width(pr_gr)+us,fix='end')

genome<-BSgenome.Hsapiens.UCSC.hg38
pr_seq<-getSeq(genome, pr_gr2)
writeFasta(pr_seq,'/Users/sonas/Desktop/motif_analysis/pr_fasta.fa')

## If you want to only look for 50 upstream of PR (excluding PR seq itself), need to update ###

# Previously, resize extended strand 50 bases upstream of the given seq (i.e. start-50 for + strand and end+50 for - strand)

# Update new coordinates to store only this upstream seq


pr_gr_df<-as.data.frame(pr_gr2)

start_coord<-ifelse(pr_gr_df$strand=='+',(pr_gr_df$start-1),(pr_gr_df$end+1))

end_coord<-ifelse(pr_gr_df$strand=='+',pr_gr_df$end,pr_gr_df$start)

pr_gr_df$start_coord<-start_coord
pr_gr_df$end_coord<-end_coord

pr_gr_df<-dplyr::rename(pr_gr_df,
                      pr_start='start',
                      pr_end='end')

pr_gr_df$start<-ifelse(pr_gr_df$strand=='+',pr_gr_df$start_coord,pr_gr_df$end_coord)

pr_gr_df$end<-ifelse(pr_gr_df$strand=='+',pr_gr_df$end_coord,pr_gr_df$start_coord)


pr_us_gr<-makeGRangesFromDataFrame(pr_gr_df,keep.extra.columns = TRUE)

genome<-BSgenome.Hsapiens.UCSC.hg38
pr_us_seq<-getSeq(genome, pr_us_gr)

writeFasta(pr_us_seq,'/Users/sonas/Desktop/motif_analysis/pr_us_fasta.fa')

