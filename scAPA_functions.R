
##### This function performs APA testing using DEXSeq and was adaped from Tinglab APA pipeline (https://github.com/2019surbhi/apa_atingLab2019/tree/master/01_polyAseq) #####

# stab: table with groups1 and group2 sample names to be tested
# peak_mat: matrix with read count per peak (row) for each sample replicate (column)
# meanmat: matrix with mean values per sample (each column is 1 sample)
# peak_ref: Final annotated peak ref table with TU (transcription unit) assiged to each peak
# 
testAPA2<-function(stab,peak_mat,meanmat,peak_ref,jtu=jtu,a=a,b=b,adjust.var=adjust.var,prs_gt=prs_gt,ncpu=20)
{

# Filter PRs based on read count within samples to be tested #

peak.expr <- rownames(meanmat)[(meanmat[,a]>=cov)|(meanmat[,b]>=cov)] 
peak.a <- rownames(meanmat)[(meanmat[,a]>=cov)]
peak.b <- rownames(meanmat)[(meanmat[,b]>=cov)]


join.a <- peak_ref
join.a <- join.a[(join.a$peak %in% peak.a)&(join.a$unique_tu==TRUE),]

join.b <- peak_ref
join.b <- join.b[(join.b$peak %in% peak.b)&(join.b$unique_tu==TRUE),]

join.keep <- peak_ref[((peak_ref$peak %in% peak.expr)&(peak_ref$unique_tu==TRUE)),]
  

stopifnot(!str_detect(join.keep$over_tus,","))
stopifnot(!str_detect(join.keep$flank_tus,","))
stopifnot(all(join.keep$unique_peak))
stopifnot(length(unique(join.keep$peak))==nrow(join.keep))
stopifnot(class(join.keep$peak)=="character")
#stopifnot(join.keep$peak %in% rownames(peak_mat))
join.keep  <- data.table(join.keep, key="tu")

join.keep <- join.keep[,list(peak=peak,
                             type=type,
                             coding=coding,
                             unique_peak=unique_peak,
                             unique_tu=unique_tu,
                             over_tus=over_tus,
                             flank_tus=flank_tus,
                             num_peak=length(peak)),
                       by="tu"]

# Remove TU with 1 peak remaining after filtering
tu_count<-join.keep$tu %>% count()
rem<-which(tu_count$freq==1)
rem_tu<-tu_count$x[rem]
rem_tu_idx<-match(rem_tu,join.keep$tu)
join.keep<-join.keep[-rem_tu_idx,]

# Now subset peak in peak_mat and meanmat
peak_mat<-peak_mat[join.keep$peak,]
meanmat<-meanmat[join.keep$peak,]
stopifnot(rownames(peak_mat)==join.keep$peak)
#rownames(mat) <- paste(join.keep$tu,join.keep$peak,sep=":")

#stab<-compdt
cd <- peak_mat[,stab$group %in% c(a,b)]
sd <- stab[stab$group %in% c(a,b),]
#sd <-colnames(peak_mat) %>% gsub('_[0-9]','',.)
#sd<-factor(sd,levels=c(a,b))
stopifnot(stab$sample==colnames(peak_mat))

# Subset peak ref based on what is found in the two comp groups
idx<-match(join.keep$peak,peak_ref$peak)
head(join.keep$peak)
peak_ref$peak[idx] %>% head()

# Create genomic range
peak<-peak_ref[idx,c('chr','start','end','strand','peak')]
peak_gr<-makeGRangesFromDataFrame(peak,keep.extra.columns = TRUE)
#fr <- peak$peak[peak$peak$peak %in% join.keep$peak]
#fr <- fr[match(join.keep$peak,fr$peak)]

stopifnot(peak_gr$peak==join.keep$peak)
rownames(cd) <- NULL

sd.use <- data.frame(group=sd$group)
sd.use$group <- factor(sd.use$group,ordered=F)

# sd.use <- data.frame(group=sd)
# sd.use$group <- factor(sd.use$group,ordered=F)


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


dxd <- DEXSeqDataSet(countData=cd, sampleData=sd.use, design=as.formula(design), featureID=join.keep$exon_num, groupID=join.keep$tu, featureRanges=pr_gr)

# Run tests
message("Estimating size factors")
dxd <- estimateSizeFactors(dxd)
message("Estimating dispersions")
dxd <- estimateDispersions(dxd)

if(is.null(adjust.var))
{
  dxd <- testForDEU(dxd)
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
counts.norm <- counts(dxd, normalized=TRUE)
counts.norm <- counts.norm[,1:nrow(sd)]
colnames(counts.norm) <- sd$sample
stopifnot(rownames(counts.norm)==join.keep$key)
rownames(counts.norm) <- join.keep$peak

# Start building nice results list
dt <- data.table(as.data.frame(dxr))
out <- with(dt,data.table(tu=groupID,peak=join.keep$peak,p=pvalue,padj=padj,dexl2fc_b_over_a=get(paste0("log2fold_",a,"_",b)),chr=genomicData.seqnames,start=genomicData.start,end=genomicData.end,strand=genomicData.strand))
out$gene_name <- red_ens$tu[match(out$tu,red_ens$tu$tu)]$name
stopifnot(paste0(out$peak)==rownames(counts.norm))

#browser()

# Make usage matrices
c2 <- data.table(cd)
c2$peak <- out$peak
c2$tu <- out$tu
cm <- melt(c2,id.vars=c("peak","tu"))
stopifnot(!is.na(cm$value))
# Within each TU and each sample, get each percent as percent of total
#m2 <- cm[,list(peak=peak,raw_count=value,use_frac=ifelse(sum(value)==0,0,value/sum(value))),by=c("tu","variable")]
cm<-as.data.table(cm)
m2 <- cm[,list(peak=peak,raw_count=value,use_frac=value/sum(value),tu_zero=sum(value)==0),by=c("tu","variable")]
stopifnot(nrow(m2[is.na(use_frac),])==sum(m2$tu_zero))
m2[is.na(use_frac),use_frac:=0]

stopifnot(!is.na(m2$use_frac))
m2$group <- stab[match(m2$variable,stab$sample),]$group
#m2$group<-stri_sub(m2$variable, 1, -3)
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

#browser()

#myraw <- cd
#colnames(myraw) <- paste0(colnames(cd),"_raw")
myuse <- use.frac
colnames(myuse) <- paste0(colnames(use.frac),"_frac")
mynorm <- counts.norm
colnames(mynorm) <- paste0(colnames(mynorm),"_norm")

nmeans_a <- rowMeans(counts.norm[,as.character(stab[stab$group==a,]$sample)])
nmeans_b <- rowMeans(counts.norm[,as.character(stab[stab$group==b,]$sample)])

# nmeans_a <- rowMeans(counts.norm[,c('u4','u4')])
# nmeans_b <- rowMeans(counts.norm[,c('b2','b2')])

dt <- data.table(mean_a=nmeans_a,mean_b=nmeans_b,b_minus_a=nmeans_b-nmeans_a,l2_b_over_a=log2((nmeans_b+0.0001)/(nmeans_a+0.0001)))
setnames(dt,paste0(colnames(dt),"_norm"))

res <- cbind(out,mynorm,dt,myuse,m3[,c("mean_a_frac","mean_b_frac","b_minus_a_frac","l2_b_over_a_frac"),with=F])

rownames(cd) <- res$peak

#browser()

write.xlsx(res,paste0(outdir,'Differential_usage.xlsx'))
write.xlsx(use.frac,paste0(outdir,'Frac_usage.xlsx'))

write.xlsx(res,paste0(outdir,'Differential_usage_w_adj.xlsx'))
write.xlsx(use.frac,paste0(outdir,'Frac_usage_w_adj.xlsx'))

# Return
res<-list(dxd=dxd,dxr=dxr,res=res,counts.raw=cd,counts.norm=counts.norm,use.frac=use.frac,a=a,b=b,stab.sub=sd)

return(res)
}
