
# Read Seurat object - using ureter paper's annotation
obj<-readRDS('/home/sonas/tingalab/Manuscripts/2021_UreterManuscript/ANALYSIS/Seurat/Uro subset/2021_04_01_ureter10_uro_PC50_res0.2/2021_04_01_ureter10_uro_PC50_res0.2_clustered.rds')

# Get metadata
meta<-obj@meta.data

# Add sample name as prefix and remove number suffix
rownames(meta)<-paste0(meta$orig.ident,'_',rownames(meta))
rownames(meta)<-stri_sub(rownames(meta), 1, -3)
rownames(meta)<-gsub('1_','1',rownames(meta)) # in case suffix is 2 digit number

bc<-rownames(meta)
meta<-cbind(bc,meta)

write.xlsx(meta,paste0(outdir,'2021_04_01_ureter10_uro_PC50_res0.2_meta.xlsx'))

## (B) Create Merged counts - need to run once per dataset ##

count_dir<-'/home/sonas/tingalab/Surbhi/PROJECTS_tinglab_drive/APA_Project/Ureter_uro/5.Seurat_workflow/counts_mat/l300iA10mm3/'

f<-list.files(count_dir,full.names = TRUE)
peak_counts<-lapply(f,read.table,sep='\t', header=TRUE,stringsAsFactors=FALSE)


# Replace . with - in count mat 
for(i in 1:length(peak_counts))
{
  colnames(peak_counts[[i]])<-gsub('.','-',colnames(peak_counts[[i]]),fixed = TRUE)
  colnames(peak_counts[[i]])<-gsub('uro_','',colnames(peak_counts[[i]]))
}

# Merge counts
merged_counts<-do.call(cbind,peak_counts)

write.table(merged_counts,paste0(outdir,'u10_peak_counts.txt'),sep='\t')

