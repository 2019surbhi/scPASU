# This code is modified from Ting lab's bulk APA pipeline to repurpose it for scRNA APA analysis
# Author: Surbhi Sona

# original script from bulkAPA pipeline: '08_pas_stats_figure.r'
# hg38 polya db: atlas.clusters.2.0.GRCh38.96.bed  (downloaded from :https://polyasite.unibas.ch/atlas#2)
# hg19 db: polyAsite_v1_hg19_clusters.bed


library(data.table)
library(goldmine)
library(handy)
library(IRanges)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Hsapiens.UCSC.hg19)


countHexes <- function(gr,us=50)
{
  hexes <- c("AATAAA","ATTAAA","AGTAAA","TATAAA","CATAAA","GATAAA","AATATA","AATACA","AATAGA","AATGAA","ACTAAA","AACAAA","TTTAAA")
  
  message("Generating flanks")
  
  # Also including the PR itself in addition to the flank here
  fl <- resize(gr,us+width(gr),fix="end")
  #fl <- flank(re,width=us)
  
                   
  message("Pulling sequence data from hg38")
  seq <- getSeq(BSgenome.Hsapiens.UCSC.hg38,names=fl)

  message("Searching hexamer occurances")
  pd <- PDict(x=hexes,max.mismatch=0)
  vc <- vcountPDict(pd, seq)
  vc <- t(vc)
  
  
  message("Selecting matches")
  ap <- apply(vc,1,function(x) which(x>=1)[1])
  hp <- hexes[ap]
  gr$hex <- hp
  gr[rowSums(vc)==0]$hex <- "NONE"
  #table(gr$hex)
  
  gr$hex <- factor(gr$hex,levels=c(hexes,"NONE"))
  
  return(gr)
}

countHexes_up <- function(gr,us=50)
{
  hexes <- c("AATAAA","ATTAAA","AGTAAA","TATAAA","CATAAA","GATAAA","AATATA","AATACA","AATAGA","AATGAA","ACTAAA","AACAAA","TTTAAA")
  
  message("Generating flanks")
  
# Save original coord
df1<-as.data.frame(gr)

# Also including the PR itself in addition to the flank here
fl <- resize(gr,us+width(gr),fix="end")
#fl <- flank(re,width=us)

# extended coord
df2<-as.data.frame(fl)

# only 50 base upstream region
df3<-df1
df3$start<-ifelse(df3$strand=='+',df2$start,df1$end)
df3$end<-ifelse(df3$strand=='+',df1$start,df2$end)

fl2<-makeGRangesFromDataFrame(df3,keep.extra.columns = TRUE)


message("Pulling sequence data from hg38")
seq2 <- getSeq(BSgenome.Hsapiens.UCSC.hg38,names=fl2)

message("Searching hexamer occurances")
pd <- PDict(x=hexes,max.mismatch=0)
vc <- vcountPDict(pd, seq2)
vc <- t(vc)


message("Selecting matches")
ap <- apply(vc,1,function(x) which(x>=1)[1])
hp <- hexes[ap]
gr$hex <- hp
gr[rowSums(vc)==0]$hex <- "NONE"
#table(gr$hex)

gr$hex <- factor(gr$hex,levels=c(hexes,"NONE"))

return(gr)
}


countHexes2 <- function(gr,us=50)
{
  hexes <- c("AATAAA","ATTAAA","AGTAAA","TATAAA","CATAAA","GATAAA","AATATA","AATACA","AATAGA","AATGAA","ACTAAA","AACAAA","TTTAAA")
  
  message("Generating flanks")
  
  # Also including the PR itself in addition to the flank here
  fl <- resize(gr,us+width(gr),fix="end")
  #fl <- flank(re,width=us)
  
  message("Pulling sequence data from hg19")
  seq <- getSeq(BSgenome.Hsapiens.UCSC.hg19,names=fl)
  
  message("Searching hexamer occurances")
  pd <- PDict(x=hexes,max.mismatch=0)
  vc <- vcountPDict(pd, seq)
  vc <- t(vc)
  
  message("Selecting matches")
  ap <- apply(vc,1,function(x) which(x>=1)[1])
  hp <- hexes[ap]
  gr$hex <- hp
  gr[rowSums(vc)==0]$hex <- "NONE"
  #table(gr$hex)
  
  gr$hex <- factor(gr$hex,levels=c(hexes,"NONE"))
  
  return(gr)
}
