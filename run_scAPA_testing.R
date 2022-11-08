source('scAPA_functions.R')

# User inputs
outdir<-'/home/sonas/APA/scAPA/'
output_file<-'u10_uro_apa.rd'
counts_file<-'u10_peak_counts.txt'

# Read merged counts
merged_counts<-read.table(paste0(outdir,counts_file),sep='\t')
colnames(merged_counts)<-gsub('.','-',colnames(merged_counts),fixed=TRUE)

# Read meta data
meta<-read.xlsx(paste0(outdir,'2021_04_01_ureter10_uro_PC50_res0.2_meta.xlsx'))

# Read peak ref
peak_ref<-read.table('/home/sonas/tingalab/Surbhi/PROJECTS_tinglab_drive/APA_Project/Ureter_uro/4.assign_tu/l300i_A10mm3_3polyA/u10_uro_final_peak_ref_l300i_A10mm3_3polya_updated.txt',header=TRUE)

# Replace peak name column with final annotation as this is more meaningful name 
col<-which(colnames(peak_ref)=='final_annotation')
peak_ref$peak<-peak_ref[,col]
peak_ref<-peak_ref[,-col]

# Remove P0 peaks

plist<-strsplit(peak_ref$pr,split=':')
tlist<-sapply(plist,'[',3)
rem<-which(tlist=='P0') # 4272 P0 peaks
peak_ref<-peak_ref[-rem,]

# Read sample and comparison tables #

comp_file<-'/home/sonas/APA/scAPA/comp_group.csv'
compdt <- fread(comp_file,header=TRUE)
comps <- data.table(a=compdt$group1,b=compdt$group2)

compnames <- paste0(comps$a,"_vs_",comps$b)

# Create stab

group<-c(comps$a,comps$b) %>% unique() %>% sort()
group<-rep(group,each=2)
sample<-c(paste0(group,'_1'),paste0(group,'_2')) %>% sort()
batch<-rep(c(1,2),(length(group)/2))
stab<-cbind(sample,group,batch) %>% as.data.frame()

stab$group<- stab$group %>% as.factor()
stab$batch<- stab$batch %>% as.factor()

cov<-10


## Run APA testing ##

dotests <- function(a,b,ncpu=4,min_peak=1)
{
	message("Testing: ",a," vs ",b)
	mycomp <- paste0(a,"_vs_",b)
	# Run without batch correction
  ta <- testApa2(stab=stab,meanmat=meanmat,peak_ref=peak_ref,peak_mat=peak_mat,jtu=jtu,a=a,b=b,adjust.var=NULL,ncpu=ncpu,min_peak=min_peak)

  # Run with batch correction
  taadj <- testApa2(stab=stab,meanmat=meanmat,peak_ref=peak_ref,peak_mat=peak_mat,jtu=jtu,a=a,b=b,adjust.var=adjust.var,ncpu=ncpu,min_peak=min_peak)
	list(ta=ta,taadj=taadj)

}

compnames <- paste0(comps$a,"_vs_",comps$b)
apa <- lapply(1:nrow(comps),function(x) dotests(a=comps[x,]$a,b=comps[x,]$b,min_peak=min_peak))
names(apa) <- compnames

# Merge down the results to one table going forward
handy::nsummary(apa)

# Put the adj info into the main res so get one ta object for each

tas <- lapply(apa,function(x){
	x$taadj$res$noadj_p <- x$ta$res$p
	x$taadj$res$noadj_padj <- x$ta$res$padj

	x$taadj$res$batchadj_p <- x$taadj$res$p
	x$taadj$res$batchadj_padj <- x$taadj$res$padj

	x$taadj$res[,p:=NULL]
	x$taadj$res[,padj:=NULL]

	return(x$taadj)
})

apa <- tas

save(apa,compress=T,file=output_file)
	      
## 2. Run t test ##
nrepa<-2
nrepb<-2

apa.tt<-testTTest(apa,nrepa,nrepb)
write.xlsx(apa.tt$res,paste0(outdir,'Differential_usage_ttest.xlsx'))

## 3. Run EdgeR ##
edger_prefix <- file.path(outdir,'edger_')
apa.alt <- testEdgeR(apa.tt)

write.xlsx(apa.alt$res,paste0(outdir,'Differential_usage_edger.xlsx'))
save(apa.alt,file=paste0(outdir,'apa.alt.rd'),compress=T)

## 4. Check significance ##
apa.sig<-callsig(apa.alt)
save(apa.sig,file=paste0(outdir,'apa.sig.rd'),compress=T)

tab<-apa.sig$res
write.xlsx(tab,paste0(outdir,'Differential_usage_w_sig.xlsx'))



