library(MURP)
library(Seurat)
library(dplyr)
library(stringi)
library(caret)

source('/home/sonas/scripts/final_scripts/my_scRNA_pipeline_functions_developing.R')

uro<-readRDS('/home/sonas/tingalab/Manuscripts/2021_UreterManuscript/ANALYSIS/Seurat/Uro subset/2021_04_01_ureter10_uro_PC50_res0.2/2021_04_01_ureter10_uro_PC50_res0.2_clustered.rds')

# Read peak obj
pobj<-readRDS('/home/sonas/tingalab/Surbhi/PROJECTS_tinglab_drive/APA_Project/scPASU/output/4b_scPASU_clustering/Seurat/peak_mat/2023_03_011_P0_500_Pn_5000_PC40res0.2/2023_03_011_P0_500_Pn_5000_PC40res0.2_clustered.rds')

# Sanity checks
identical(rownames(uro@meta.data),colnames(uro))
identical(rownames(pobj@meta.data),colnames(pobj))

# Correct barcode in peak obj
new_bc<-gsub('\\.','-',rownames(pobj@meta.data))
new_bc<-gsub('_uro','',new_bc)
pobj<-RenameCells(pobj,new.names = new_bc)

# correct seurat obj barcode
new_bc_s<-paste0(uro$orig.ident,'_',rownames(uro@meta.data))
new_bc_s<-stri_sub(new_bc_s, 1, -3)
new_bc_s<-gsub('1_','1',new_bc_s) # in case there is more than 9 samples
uro<-RenameCells(uro,new.names = new_bc_s)


# Transfer Seurat clusters
idx<-match(rownames(uro@meta.data),rownames(pobj@meta.data))
pobj@meta.data[idx,"old_clusters"]<-uro@meta.data$seurat_clusters


clus_obj <- SplitObject(pobj, split='old_clusters') # seurat_clusters column in metadata stores my cluster identity

clus_obj<-clus_obj[order(names(clus_obj))]

clusters<-c('i0','i1','b2','b3','u4','i5','b6','u7')

cols<-c('#FABED4','#BFEF45','#FF0000','#FFBF00','#469990','#DF00FF','#800000','#000075')

#### Data partitioning ###

outdir<-'/home/sonas/thesis_figures/ch4/data_partition/umaps/'
rep_dir<-'/home/sonas/thesis_figures/ch4/data_partition/rep_dir/'

for(i in 1:length(clusters))
{
# extract normalized counts
clus <- clus_obj[[i]]@assays$RNA@data %>% as.matrix()
input <- t(clus)

# 2. running murp
m= MURP(Data = input, max_murp = 2, omega = 1/6, 
        cores = 30, iter=5,fast = T,seed = 01234567)

# Tried k=3  with 4 iter for the two umbrella clusters
m= MURP(Data = input, max_murp = 3, omega = 1/6, 
        cores = 30, iter = 4,fast = T,seed = 01234567)

r1_bc<-names(m$Recommended_K_cl$cluster)[m$Recommended_K_cl$cluster==1]

r2_bc<-names(m$Recommended_K_cl$cluster)[m$Recommended_K_cl$cluster==2]

# Create replicates

#rep1<-subset(pobj,cells = r1_bc)
#rep2<-subset(pobj,cells = r2_bc)

write.table(as.data.frame(r1_bc),paste0(rep_dir,clusters[i],'_1.tsv'),sep='\t',row.names = FALSE,col.names = FALSE,quote=FALSE)

write.table(as.data.frame(r2_bc),paste0(rep_dir,clusters[i],'_2.tsv'),sep='\t',row.names = FALSE,col.names = FALSE,quote=FALSE)



png(paste0(outdir,'cluster_',clusters[i],'.png'),width = 10,height = 8,units = "in",res = 600)
DimPlot(pobj,group.by = 'old_clusters',cells.highlight =  c(r1_bc,r2_bc),
        cols.highlight = cols[i],pt.size = 0.5)
dev.off()

png(paste0(outdir,'cluster_',clusters[i],'_rep1.png'),width = 10,height = 8,units = "in",res = 600)
DimPlot(pobj,group.by = 'old_clusters',cells.highlight =  r1_bc,
        cols.highlight =cols[i],pt.size = 0.5)
dev.off()

png(paste0(outdir,'cluster_',clusters[i],'_rep2.png'),width = 10,height = 8,units = "in",res = 600)
DimPlot(pobj,group.by = 'old_clusters',cells.highlight =  r2_bc,
        cols.highlight = cols[i],pt.size = 0.5)
dev.off()



### URO UMAP
png(paste0(outdir,'uro_cluster_',clusters[i],'.png'),width = 10,height = 8,units = "in",res = 600)
DimPlot(uro,group.by = 'seurat_clusters',cells.highlight =  c(r1_bc,r2_bc),
        cols.highlight = cols[i],pt.size = 0.5)
dev.off()

png(paste0(outdir,'uro_cluster_',clusters[i],'_rep1.png'),width = 10,height = 8,units = "in",res = 600)
DimPlot(uro,group.by = 'seurat_clusters',cells.highlight =  r1_bc,
        cols.highlight =cols[i],pt.size = 0.5)
dev.off()

png(paste0(outdir,'uro_cluster_',clusters[i],'_rep2.png'),width = 10,height = 8,units = "in",res = 600)
DimPlot(uro,group.by = 'seurat_clusters',cells.highlight =  r2_bc,
        cols.highlight = cols[i],pt.size = 0.5)
dev.off()


}

### Had issues with partitioning the two umbrella clusters.





