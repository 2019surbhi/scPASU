source('scAPA_functions.R')

# User inputs
outdir<-'/home/sonas/APA/scAPA/'
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



