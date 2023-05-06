library(data.table)
library(goldmine)
library(handy)
library(IRanges)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Hsapiens.UCSC.hg19)


polyadb_fpath_hg38<-'/home/sonas/tingalab/Surbhi/PROJECTS_tinglab_drive/APA_Project/scPASU/input/atlas.clusters.2.0.GRCh38.96.bed'

f<-'/home/sonas/tingalab/Surbhi/PROJECTS_tinglab_drive/APA_Project/scPASU/output/2_PeakRef/2e_TSS_filtered/u10_uro_l300iA10mm3_e100_updated_tss_filt.txt'

out_dir<-'/home/sonas/APA/output/hexamer_plot/output/'
tmp_dir<-'/home/sonas/APA/output/hexamer_plot/tmp/'

file<-read.table(f,header=TRUE)

if(class(file)=='data.frame')
{
  gr<-if(class(file)=='data.frame')
  {
    gr<-makeGRangesFromDataFrame(file)
    
  }else{
    gr<-file
  }(file)
  
}else{
  gr<-file
}


### inferred polyA site ###
file2<-file
colnames(file2)<-gsub('\\bchr\\b','peak_chr',colnames(file2))
colnames(file2)<-gsub('\\bstart\\b','peak_start',colnames(file2))

colnames(file2)<-gsub('\\bend\\b','peak_end',colnames(file2))
colnames(file2)<-gsub('\\bstrand\\b','peak_strand',colnames(file2))

colnames(file2)<-gsub('\\bpr_chr\\b','chr',colnames(file2))
colnames(file2)<-gsub('\\bpr_start\\b','start',colnames(file2))
colnames(file2)<-gsub('\\bpr_end\\b','end',colnames(file2))
colnames(file2)<-gsub('\\bpr_strand\\b','strand',colnames(file2))

gr2<-makeGRangesFromDataFrame(file2)


### pdb 
chr<-paste0('chr',1:22)
chr<-c(chr,'chrX','chrY')

pdb_hg38 <- fread(polyadb_fpath_hg38)
pdb_hg38$V1<-paste0('chr',pdb_hg38$V1)
pdb_hg38<-pdb_hg38[pdb_hg38$V1 %in% chr,]

# Create GR range
pa_hg38 <- with(pdb_hg38,GRanges(V1,IRanges(V2,V3),strand=V6))
pa_hg38 <- pa_hg38[seqnames(pa_hg38) %in% handy::chrs()]
pa_hg38 <- unique(pa_hg38)


# # Extend by 15 bases on either side
pa_hg38<-resize(pa_hg38,width = width(pa_hg38)+15,fix = 'start')
pa_hg38<-resize(pa_hg38,width = width(pa_hg38)+15,fix = 'end')


## Look for overlap ##

# pr contains 58,191 PR
# polyA db contains 568,608 pA

library(ChIPpeakAnno)
library(Vennerable)
library(ggplot2)

res<-makeVennDiagram(list(gr2, pa_hg38), NameOfPeaks=c("PR", "polyASite"), minoverlap = 1,ignore.strand = FALSE,
                     scaled=FALSE, euler.d=FALSE,
                     fill=c("#009E73", "#F0E442"), # circle fill color
                     col=c("#D55E00", "#0072B2"), #circle border color
                     cat.col=c("#D55E00", "#0072B2"))





venn_cnt2venn <- function(venn_cnt){
  n <- which(colnames(venn_cnt)=="Counts") - 1
  SetNames=colnames(venn_cnt)[1:n]
  Weight=venn_cnt[,"Counts"]
  names(Weight) <- apply(venn_cnt[,1:n],1, paste, collapse="")
  Venn(SetNames=SetNames,Weight = Weight )
}

v <- venn_cnt2venn(res$vennCounts)

out<-'/home/sonas/APA/output/hexamer_plot/output/'
png(paste0(out,'overlap_pr_polyASite.png'))
plot(v)
dev.off()





res2<-makeVennDiagram(list(gr, pa_hg38_clus), NameOfPeaks=c("peak", "polyASite"), minoverlap = 1,ignore.strand = FALSE,
                      scaled=FALSE, euler.d=FALSE,
                      fill=c("#009E73", "#F0E442"), # circle fill color
                      col=c("#D55E00", "#0072B2"), #circle border color
                      cat.col=c("#D55E00", "#0072B2"))



v2 <- venn_cnt2venn(res2$vennCounts)

png(paste0(out,'overlap_peak_polyASite.png'))
plot(v2)
dev.off()
