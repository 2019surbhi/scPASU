library(dplyr)
library(stringi)
library(stringr)
library(readxl)
library(openxlsx)
library(data.table)
library(matrixStats)
library(gtools)

library(ggplot2)
library(reshape)

library(Seurat)
library(Rsamtools)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Hsapiens.UCSC.hg38)
library(edgeR)
library(DEXSeq)


##### This function performs APA testing using DEXSeq and was adaped from Tinglab APA pipeline (https://github.com/2019surbhi/apa_atingLab2019/tree/master/01_polyAseq) #####

# stab: table with groups1 and group2 sample names to be tested
# peak_mat: matrix with read count per peak (row) for each sample replicate (column)
# meanmat: matrix with mean values per sample (each column is 1 sample)
# peak_ref: Final annotated peak ref table with TU (transcription unit) assiged to each peak
# jtu:
# a: 1st sample to be tested
# b: 2nd sample name to be tested
# adjust.var: variable to adjust for in statistical testing. Default is 'batch'
# min_peak: threshold for retaining TU with peaks greater than this number. Deafult =1 
#

testAPA2<-function(stab,peak_mat,meanmat,peak_ref,a=a,b=b,adjust.var=NULL,min_peak,ncpu=20)
{
  
  message("Creating inputs for DEXSeq")
  
  # Filter peaks based on read count within samples to be tested #
  peak.expr <- rownames(meanmat)[(meanmat[,a]>=cov)|(meanmat[,b]>=cov)]
  
  # Retain only those peaks in this new peak ref that made the read cutoff (default 10) as well as assigned to unique TU (Transcription Unit) #
  join.keep <- peak_ref[((peak_ref$peak %in% peak.expr)&(peak_ref$unique_tu==TRUE)),]
  
  # Sanity checks #
  stopifnot(!str_detect(join.keep$over_tus,","))
  stopifnot(!str_detect(join.keep$flank_tus,","))
  stopifnot(all(join.keep$unique_peak)) 
  stopifnot(length(unique(join.keep$peak))==nrow(join.keep))
  stopifnot(class(join.keep$peak)=="character")
  stopifnot(join.keep$peak %in% rownames(peak_mat))
  
  # Create data table for this new peak ref #
  join.keep  <- data.table(join.keep, key="tu")
  
  # Add relevant columns by TU #
  join.keep <- join.keep[,list(peak=peak,
                               type=type,
                               coding=coding,
                               unique_peak=unique_peak,
                               unique_tu=unique_tu,
                               over_tus=over_tus,
                               flank_tus=flank_tus,
                               num_peak=length(peak)),
                         by="tu"]
  
  # Remove TU with 1 peak remaining after applying all filters #
  join.keep <- join.keep[num_peak > min_peak,] 
  
  #tu_count<-join.keep$tu %>% count()
  #rem<-which(tu_count$freq==1)
  #rem_tu<-tu_count$x[rem]
  #rem_tu_idx<-match(rem_tu,join.keep$tu)
  #join.keep<-join.keep[-rem_tu_idx,]
  
  # Now subset peak in peak_mat and meanmat to include only the peaks from this new peak ref #
  peak_mat<-peak_mat[join.keep$peak,]
  meanmat<-meanmat[join.keep$peak,]
  
  stopifnot(rownames(peak_mat)==join.keep$peak)
  
  sd <- stab[stab$group %in% c(a,b),] %>% arrange(.,group)
  
  # Sort by peak mat column
  stopifnot(sd$sample %in% colnames(peak_mat))
  sd<-sd[match(sd$sample,colnames(peak_mat)),]
  cd <- peak_mat
  
  #cd <- peak_mat[,stab$group %in% c(a,b)]
  # cols<-colnames(cd) %>% sort()
  # cd<-cd[,cols]
  
  
  # Subset peak ref based on what is found in the two comp groups
  idx<-match(join.keep$peak,peak_ref$peak)
  
  # Create genomic range
  peak<-peak_ref[idx,c('chr','start','end','strand','peak')]
  peak_gr<-makeGRangesFromDataFrame(peak,keep.extra.columns = TRUE)
  
  stopifnot(peak_gr$peak==join.keep$peak)
  sd$group <- factor(sd$group,levels=c(a,b))
  rownames(cd) <- NULL
  
  sd.use <- data.frame(group=sd$group)
  sd.use$group <- factor(sd.use$group,ordered=F)
  
  
  message("Building DEXSeqDataSet object")
  if(is.null(adjust.var))
  {
    design <- "~ sample + exon + group:exon"
  }else
  {
    design <- paste0("~ sample + exon + adjust:exon + group:exon")
    sd.use$adjust <- factor(with(sd,get(adjust.var)))
  }
  message("Design Formula: ",design)
  
  exon_num <- join.keep[,list(paste0("E",1:length(peak)),peak=peak),by="tu"]
  stopifnot(exon_num$tu==join.keep$tu)
  stopifnot(exon_num$peak==join.keep$peak)
  join.keep$exon_num <- exon_num$V1
  join.keep$key <- paste(join.keep$tu,join.keep$exon_num,sep=":")
  
  
  dxd <- DEXSeqDataSet(countData=cd, sampleData=sd.use, design=as.formula(design), featureID=join.keep$exon_num, groupID=join.keep$tu, featureRanges=peak_gr)
  
  # Run tests
  message("Estimating size factors")
  dxd <- DEXSeq::estimateSizeFactors(dxd)
  message("Estimating dispersions")
  dxd <- DEXSeq::estimateDispersions(dxd)
  
  if(is.null(adjust.var))
  {
    message("Tesing w/o adjustment")
    dxd<-testForDEU(dxd)
  }else
  {
    message("Tesing w/ adjustment for variable: ", adjust.var)
    dxd <- testForDEU(dxd,fullModel=as.formula("~ sample + exon + adjust:exon + group:exon"),reducedModel=as.formula("~ sample + exon + adjust:exon"))
  }
  
  
  dxd <- estimateExonFoldChanges(dxd,fitExpToVar="group",BPPARAM=MulticoreParam(workers=ncpu))
  
  # Pull APA results table
  message("Extracting DEXSeqResults")
  dxr <- DEXSeqResults(dxd)
  
  # Fetch normalized counts
  counts.norm <- DEXSeq::counts(dxd, normalized=TRUE)
  counts.norm <- counts.norm[,1:nrow(sd)]
  colnames(counts.norm) <- sd$sample
  stopifnot(rownames(counts.norm)==join.keep$key)
  rownames(counts.norm) <- join.keep$peak
  
  # Start building nice results list
  dt <- data.table(as.data.frame(dxr))
  out <- with(dt,data.table(tu=groupID,peak=join.keep$peak,p=pvalue,padj=padj,dexl2fc_b_over_a=get(paste0("log2fold_",b,"_",a)),chr=genomicData.seqnames,start=genomicData.start,end=genomicData.end,strand=genomicData.strand))
  
  # Extract gene name
  genes<-strsplit(out$peak,split=':') %>% sapply(.,'[[',2) 
  out$gene_name<-genes
  stopifnot(out$peak==rownames(counts.norm))
  
  # Make usage matrices
  c2 <- data.table(cd)
  c2$peak <- out$peak
  c2$tu <- out$tu
  cm <- melt(c2,id.vars=c("peak","tu"))
  stopifnot(!is.na(cm$value))
  
  # Within each TU and each sample, get each percent as percent of total
  cm<-as.data.table(cm)
  m2 <- cm[,list(peak=peak,raw_count=value,use_frac=value/sum(value),tu_zero=sum(value)==0),by=c("tu","variable")]
  stopifnot(nrow(m2[is.na(use_frac),])==sum(m2$tu_zero))
  
  m2[is.na(use_frac),use_frac:=0]
  stopifnot(!is.na(m2$use_frac))
  
  m2$group <- stab[match(m2$variable,stab$sample),]$group
  stopifnot(!is.na(m2$group))
  
  # Groupwise usage means
  m3 <- m2[,list(tu=tu[1],mean_a=mean(use_frac[group==a]),mean_b=mean(use_frac[group==b])),by=c("peak")]
  m3$b_minus_a <- m3$mean_b-m3$mean_a
  m3$l2_b_over_a <- log2((m3$mean_b+0.0001)/(m3$mean_a+0.0001))
  stopifnot(out$tu==m3$tu)
  stopifnot(out$peak==m3$peak)
  setnames(m3,paste0(colnames(m3),"_frac"))
  
  cas <- cast(m2,formula="peak+tu~variable",value="use_frac")
  caskey <- paste(cas$tu,cas$peak,sep=":")
  outkey <- paste(out$tu,out$peak,sep=":")
  m <- match(outkey,caskey)
  stopifnot(!is.na(m))
  cas <- cas[m,]
  use.frac <- cas[,c(-1,-2)]
  stopifnot(nrow(use.frac)==nrow(out))
  rownames(use.frac) <- join.keep$peak
  
  
  myuse <- use.frac
  colnames(myuse) <- paste0(colnames(use.frac),"_frac")
  mynorm <- counts.norm
  colnames(mynorm) <- paste0(colnames(mynorm),"_norm")
  
  nmeans_a <- rowMeans(counts.norm[,as.character(stab[stab$group==a,]$sample)])
  nmeans_b <- rowMeans(counts.norm[,as.character(stab[stab$group==b,]$sample)])
  
  dt <- data.table(mean_a=nmeans_a,mean_b=nmeans_b,b_minus_a=nmeans_b-nmeans_a,l2_b_over_a=log2((nmeans_b+0.0001)/(nmeans_a+0.0001)))
  setnames(dt,paste0(colnames(dt),"_norm"))
  
  res <- cbind(out,mynorm,dt,myuse,m3[,c("mean_a_frac","mean_b_frac","b_minus_a_frac","l2_b_over_a_frac"),with=F])
  
  rownames(cd) <- res$peak
  
  #write.xlsx(res,paste0(outdir,'Differential_usage.xlsx'))
  #write.xlsx(use.frac,paste0(outdir,'Frac_usage.xlsx'))
  
  #write.xlsx(res,paste0(outdir,'Differential_usage_w_adj.xlsx'))
  #write.xlsx(use.frac,paste0(outdir,'Frac_usage_w_adj.xlsx'))
  
  # Return
  ret<-list(dxd=dxd,dxr=dxr,res=res,counts.raw=cd,counts.norm=counts.norm,use.frac=use.frac,a=a,b=b,stab.sub=sd)
  
  return(ret)
}


##### This function is borrowed from Ting lab APA pipeline to perfor differential testing using t-test #####

# myta: apa obj with non adjusted and adjusted p val from DEXseq testing
# nrepa: n rep in group a (or group 1)
# nrepb: n rep in group b (or group 2)

testTTest <- function(myta,nrepa,nrepb)
{
  dat <- myta$use.frac
  # Found cases where both groups have zero variance and t-test failed on them, so have to flag this and pass 0 p-value
  t.p <- lapply(1:nrow(dat),function(x){
    #message(x);
    #avec <- dat[x,][,1:4]
    #bvec <- dat[x,][,5:8]
    
    # We have only 2 rep
    avec <- dat[x,][,1:nrepa]
    bvec <- dat[x,][,(nrepa+1):nrepb]
    if((sum(avec==max(avec))==nrepa)&(sum(bvec==max(bvec))==nrepb))
    {
      return(0)
    } else
    {
      td <- t.test(x=avec,y=bvec)
      return(td$p.value)
    }
  })
  t.p <- do.call(c,t.p)
  t.padj <- p.adjust(t.p,method="fdr")
  myta$res$ttest_p <- t.p
  myta$res$ttest_padj <- t.padj
  return(myta)
}

##### This function is borrowed from tinglab APA pipeline to perform differential testing using EdgeR #####

# myta: apa object with DEXSeq and t-test results
  
# Run function
testEdgeR <- function(myta)
{
  mycomp<-paste0(unique(myta$stab.sub$group)[1],"_vs_",unique(myta$stab.sub$group)[2])
  message("Testing: ",mycomp)
  
  res <- myta$res
  ann <- data.frame(tu=res$tu,pr=res$peak)
  g <- as.character(myta$stab.sub$group)
  adj <- factor(myta$stab.sub$batch)
  cnt <- myta$counts.raw
  y.all <- DGEList(counts=cnt,genes=ann,group=g)
  y <- calcNormFactors(y.all)
  design <- model.matrix(~ g + adj)
  y <- estimateDisp(y, design, robust=TRUE)
  fit <- glmQLFit(y, design, robust=TRUE)
  
  #pdf(file=paste0("output/edger_",mycomp,".pdf"))
  pdf(file=paste0(edger_prefix,mycomp,".pdf"))
  plotBCV(y)
  plotQLDisp(fit)
  plotMDS(y)
  dev.off()
  
  qlf <- glmQLFTest(fit, coef=2)
  sp <- diffSpliceDGE(fit, coef=2, geneid="tu", exonid="pr")
  
  myta$edger <- sp
  
  # Not all PR in sp
  idx<-match(names(sp$exon.p.value),myta$res$peak)
  myta$res$edger_exon_p<-rep(NA,nrow(myta$res))
  myta$res$edger_exon_p[idx] <- sp$exon.p.value
  myta$res$edger_exon_padj[idx] <- p.adjust(sp$exon.p.value,method="fdr")
  
  ftest <- sp$gene.p.value
  ftest.adj <- ftest
  ftest.adj[,1] <- p.adjust(ftest[,1],method="fdr")
  simes <- sp$gene.Simes.p.value
  
  myta$res$edger_ftest_p <- ftest[match(res$tu,rownames(ftest))]
  myta$res$edger_ftest_padj <- ftest.adj[match(res$tu,rownames(ftest.adj))]
  
  names(simes) <- res[match(names(simes),res$pr),]$tu
  simes.adj <- p.adjust(simes,method="fdr")
  
  myta$res$edger_simes_p <- simes[match(res$tu,names(simes))]
  myta$res$edger_simes_padj <- simes.adj[match(res$tu,names(simes.adj))]
  
  
  return(myta)
}


##### Call significance #####

callsig <- function(myta)
{
  # P-value and fraction delta
  myta$res$batchadj_delta_sig <- ((abs(myta$res$b_minus_a_frac)>=0.1)&(myta$res$batchadj_padj<0.0001))
  
  # P-value and fold change
  myta$res$batchadj_perfc_sig <- with(myta$res,(batchadj_padj<0.0001)&(abs(l2_b_over_a_frac)>log2(1.5))&((mean_a_frac>=0.05)|(mean_b_frac>=0.05)))
  
  # P-value and fold change and delta fraction (intersect set we want to use)
  myta$res$int_sig <- myta$res$batchadj_delta_sig & myta$res$batchadj_perfc_sig 
  
  return(myta)
}


##### This function filters APA table to inlcude TU entires for which at least 1 peak tested significant by both FC and pvalue #####

# tab: Final APA table 

filter_for_sig<-function(tab)
{
  filt<-tab[which(tab$both_sig==TRUE),]
  idx<-which(tab$tu %in% filt$tu)
  filt_tab<-tab[idx,]
  return(filt_tab)
  
}

##### This function creates subset of peak matrix by cluster cell barcode #####

# clus_name: name of cluster to create replicates for
# merged counts: read counts for peaks per cell (merged for all samples)
# meta: metadata dataframe

subset_peaks_by_clus<-function(clus,merged_counts,meta)
{
  ## extract cluster cell barcodes ##
  clus_idx<-which(meta$seurat_clusters==clus)
  clus_bc<-meta$bc[clus_idx]
  
  ## Subset peak counts data ##
  select<-match(clus_bc,colnames(merged_counts))
  
  # If an entire cell is lost due to various filer, remove it
  rem<-which(is.na(select)==TRUE)
  if(length(rem)!=0)
  {select<-select[-rem]}
  
  clus_peak<-merged_counts[,select]
  return(clus_peak)
}


##### This function creates replicates for each cluster by sampling barcode #####
# clus_name: name of cluster to create replicates for
# peak_mat: peak x cell matrix corresponding to that cluster
# nrep: number of replicates to create (with replacement). Default is 2
# p: proportion of total cells to sample. Default is 0.7
# seed: seed value. Default is 123

create_replicates<-function(clus_name,peak_mat,nrep=2,p=0.7,seed=123)
{
  
  bc<-colnames(peak_mat)
  ncell<-length(bc)
  
  rep<-list()
  
  set.seed(seed)
  for(i in 1:nrep)
  {
    subset_bc<-sample(bc,(p*ncell),replace = TRUE)
    subset_mat<-peak_mat[,subset_bc]
    rep[[i]]<-rowSums(subset_mat) %>% as.data.frame()
    names(rep)[[i]]<-paste0(clus_name,'_',i)
  }
  
  return(rep)
}


##### This function prepares the final APA table #####

# tab: datatable with APA results (DEXSeq,t-test, edgeR)
# g1: group1 sample tag
# g2: group2 sample tag

get_final_APA_tab<-function(tab,g1,g2)
{
setnames(tab,"chr","peak_chr")
setnames(tab,"start","peak_start")
setnames(tab,"end","peak_end")
setnames(tab,"strand","peak_strand")
setnames(tab,"mean_a_norm",sprintf("mean_%s_norm",g1))
setnames(tab,"mean_b_norm",sprintf("mean_%s_norm",g2))

setnames(tab,"b_minus_a_norm",sprintf("%s_minus_%s_norm",g2,g1))
setnames(tab,"l2_b_over_a_norm",sprintf("l2_%s_over_%s_norm",g2,g1))

setnames(tab,"mean_a_frac",sprintf("mean_%s_frac",g1))
setnames(tab,"mean_b_frac",sprintf("mean_%s_frac",g2))

setnames(tab,"b_minus_a_frac",sprintf("%s_minus_%s_frac",g2,g1))
setnames(tab,"l2_b_over_a_frac",sprintf("l2_%s_over_%s_frac",g2,g1))

setnames(tab,"batchadj_padj","dex_padj")
setnames(tab,"batchadj_delta_sig","delta_sig")
setnames(tab,"batchadj_perfc_sig","fracfc_sig")
setnames(tab,"int_sig","both_sig")

# Remove unwanted cols
tab[,c('dexl2fc_b_over_a','noadj_p','noadj_padj','batchadj_p','ttest_p','ttest_padj','edger_exon_p','edger_exon_padj','edger_ftest_p','edger_ftest_padj','edger_simes_p','edger_simes_padj'):=NULL]

# Order

cols<-colnames(tab)

attrib_cols<-c('gene_name',
               'tu',
               'peak',
               'peak_chr',
               'peak_start',
               'peak_end',
               'peak_strand')

stat_cols<-c('dex_padj',
             'delta_sig',
             'fracfc_sig',
             'both_sig')

rem_cols<-cols[!(cols %in% c(attrib_cols,stat_cols))]
  
  
setcolorder(tab,c(attrib_cols,rem_cols,stat_cols))

return(tab)

}




##### This function creates replicates of two given cluster and performs APA testing using DEXSeq wrapped in bulk APA analysis strategy published by our lab #####

# merged_counts: merged peak x cell counts for all samples
# meta: meta data dataframe
# c1 & c2: seurat clusters identity for 1st and 2nd clusters
# c1_nm & c2_nm: cluster names (this should match their names in sample sheet)

create_apa_inputs<-function(merged_counts,meta,c1,c2,c1_nm,c2_nm,out,plot_corr=FALSE)
{
  
  # Create subset peak matrix for each cluster
  c1_peaks<-subset_peaks_by_clus(c1,merged_counts,meta)
  c2_peaks<-subset_peaks_by_clus(c2,merged_counts,meta)
  
  # Create replicates # - adapt for more than 2 replicates?
  c1_rep_lst<-create_replicates(clus_name=c1_nm,peak_mat=c1_peaks,nrep=2,p=0.7,seed=123)
  c1_mat<-do.call(cbind,c1_rep_lst)
  colnames(c1_mat)<-names(c1_rep_lst)
  
  c2_rep_lst<-create_replicates(clus_name=c2_nm,peak_mat=c2_peaks,nrep=2,p=0.7,seed=123)
  c2_mat<-do.call(cbind,c2_rep_lst)
  colnames(c2_mat)<-names(c2_rep_lst)
  
  
  # Merge #
  peak_mat<-cbind(c2_mat,c1_mat) # same as mat
  
  # Correlation
  cor_mat<-cor(peak_mat)
  
  if(plot_corr==TRUE)
  {
    
    png(paste0(out,c1_nm,'_',c2_nm,'_correlation.png'),
        width = 8,height=8, units="in",res = 300)
    
    corrplot.mixed(corr = cor_mat,upper = 'pie',
                          lower='number',order='hclust',
                          tl.pos = "lt", tl.col = "black",
                          tl.offset=1, tl.srt = 0)
    
    dev.off()
  }  
  # Now get mean matrix by calculating mean for each sample #
  
  c1_mean<-apply(X =c1_mat,1,mean)
  c2_mean<-apply(X=c2_mat,1,mean)
  
  meanmat<-cbind(c2_mean,c1_mean) 
  colnames(meanmat)<-c(c2_nm,c1_nm)
  
  return(list(peak_mat,meanmat))
  
}


