#!/usr/bin/bash

#SBATCH --mail-type=END
#SBATCH --mail-user=sonas@ccf.org
#SBATCH --job-name=seurat_for_cell_mat
#SBATCH -n 35 -p xtreme
#SBATCH --mem-per-cpu=7900


#SBATCH -o seurat_for_cell_mat.out
#SBATCH -e seurat_for_cell_mat.err

module load R/4.2.0


apa_file='/home/sonas/beegfs/APA/scAPA/ureter10/misc/apa_genes.txt'

# u7 filtered
apa_file='/home/sonas/beegfs/APA/scAPA/ureter10/misc/apa_genes_filt.txt'

apa_file_top10='/home/sonas/beegfs/APA/scAPA/ureter10/misc/apa_genes_top10.txt'
apa_file_top20='/home/sonas/beegfs/APA/scAPA/ureter10/misc/apa_genes_top20.txt'
apa_file_top30='/home/sonas/beegfs/APA/scAPA/ureter10/misc/apa_genes_top30.txt'


#Rscript /home/sonas/beegfs/APA/scAPA/ureter10/misc/my_scRNA_pipeline_for_scpolya.R -i /home/sonas/beegfs/APA/scAPA/ureter10/rerun/peak_mat/l300iA10mm3_extn100_tss_filt/counts_mat_dir/ -o /home/sonas/beegfs/APA/scAPA/ureter10/rerun/peak_mat/l300iA10mm3_extn100_tss_filt/seurat_results/2022_12_15_u10_peak_mat_extn100_tss_filt_P0_200_Pn2000/ -f 2022_12_15_u10_peak_mat_extn100_tss_filt_P0_200_Pn2000 -C cca-mnn -z -c 40 -m 250 -v -y 200,2000 


# 200 P0 and rest APA genes (filtered for FC>1.5)
#Rscript /home/sonas/beegfs/APA/scAPA/ureter10/misc/my_scRNA_pipeline_for_scpolya.R -i /home/sonas/beegfs/APA/scAPA/ureter10/rerun/peak_mat/l300iA10mm3_extn100_tss_filt/counts_mat_dir/ -o /home/sonas/beegfs/APA/scAPA/ureter10/rerun/peak_mat/l300iA10mm3_extn100_tss_filt/seurat_results/2023_01_10_u10_peak_mat_extn100_tss_filt_P0_200_Pn_APA/ -f 2023_01_10_u10_peak_mat_extn100_tss_filt_P0_200_Pn_APA -C cca-mnn -z -c 40 -m 250 -v -y 200,2000 -Y ${apa_file} -u clustree,sil


# 200 P0 and rest APA genes (filtered for FC>1.5) and all u7 comaprisons excluded

#Rscript /home/sonas/beegfs/APA/scAPA/ureter10/misc/my_scRNA_pipeline_for_scpolya.R -i /home/sonas/beegfs/APA/scAPA/ureter10/rerun/peak_mat/l300iA10mm3_extn100_tss_filt/counts_mat_dir/ -o /home/sonas/beegfs/APA/scAPA/ureter10/rerun/peak_mat/l300iA10mm3_extn100_tss_filt/seurat_results/2023_01_12_u10_peak_mat_extn100_tss_filt_P0_200_Pn_APA/ -f 2023_01_12_u10_peak_mat_extn100_tss_filt_P0_200_Pn_APA -C cca-mnn -z -c 40 -m 250 -v -y 200,2000 -Y ${apa_file} -u clustree,sil


# optimal parameters for above
Rscript /home/sonas/beegfs/APA/scAPA/ureter10/misc/my_scRNA_pipeline_for_scpolya.R -i /home/sonas/beegfs/APA/scAPA/ureter10/rerun/peak_mat/l300iA10mm3_extn100_tss_filt/counts_mat_dir/ -o /home/sonas/beegfs/APA/scAPA/ureter10/rerun/peak_mat/l300iA10mm3_extn100_tss_filt/seurat_results/2023_02_13_u10_peak_mat_extn100_tss_filt_P0_200_Pn_APA/ -f 2023_02_13_u10_peak_mat_extn100_tss_filt_P0_200_Pn_APA -C cca-mnn -z -c 40 -m 250 -v -y 200,2000 -Y ${apa_file} -d 1:35 -e 0.5

# 200 P0 and top10% APA genes

#Rscript /home/sonas/beegfs/APA/scAPA/ureter10/misc/my_scRNA_pipeline_for_scpolya.R -i /home/sonas/beegfs/APA/scAPA/ureter10/rerun/peak_mat/l300iA10mm3_extn100_tss_filt/counts_mat_dir/ -o /home/sonas/beegfs/APA/scAPA/ureter10/rerun/peak_mat/l300iA10mm3_extn100_tss_filt/seurat_results/2023_01_12_u10_peak_mat_extn100_tss_filt_P0_200_apa_top10/ -f 2023_01_12_u10_peak_mat_extn100_tss_filt_P0_200_apa_top10 -C cca-mnn -z -c 40 -m 250 -v -y 200,2000 -Y ${apa_file_top10} -u clustree,sil


# 200 P0 and top 20% APA genes
#Rscript /home/sonas/beegfs/APA/scAPA/ureter10/misc/my_scRNA_pipeline_for_scpolya.R -i /home/sonas/beegfs/APA/scAPA/ureter10/rerun/peak_mat/l300iA10mm3_extn100_tss_filt/counts_mat_dir/ -o /home/sonas/beegfs/APA/scAPA/ureter10/rerun/peak_mat/l300iA10mm3_extn100_tss_filt/seurat_results/2023_01_12_u10_peak_mat_extn100_tss_filt_P0_200_apa_top20/ -f 2023_01_12_u10_peak_mat_extn100_tss_filt_P0_200_apa_top20 -C cca-mnn -z -c 40 -m 250 -v -y 200,2000 -Y ${apa_file_top20} -u clustree,sil

# 200 P0 and top 30% APA genes
#Rscript /home/sonas/beegfs/APA/scAPA/ureter10/misc/my_scRNA_pipeline_for_scpolya.R -i /home/sonas/beegfs/APA/scAPA/ureter10/rerun/peak_mat/l300iA10mm3_extn100_tss_filt/counts_mat_dir/ -o /home/sonas/beegfs/APA/scAPA/ureter10/rerun/peak_mat/l300iA10mm3_extn100_tss_filt/seurat_results/2023_01_12_u10_peak_mat_extn100_tss_filt_P0_200_apa_top30/ -f 2023_01_12_u10_peak_mat_extn100_tss_filt_P0_200_apa_top30 -C cca-mnn -z -c 40 -m 250 -v -y 200,2000 -Y ${apa_file_top30} -u clustree,sil



# 200 P0 and 4000 pn genes
#Rscript /home/sonas/beegfs/APA/scAPA/ureter10/misc/my_scRNA_pipeline_for_scpolya.R -i /home/sonas/beegfs/APA/scAPA/ureter10/rerun/peak_mat/l300iA10mm3_extn100_tss_filt/counts_mat_dir/ -o /home/sonas/beegfs/APA/scAPA/ureter10/rerun/peak_mat/l300iA10mm3_extn100_tss_filt/seurat_results/2023_01_11_u10_peak_mat_extn100_tss_filt_P0_200_Pn4000/ -f 2023_01_11_u10_peak_mat_extn100_tss_filt_P0_200_Pn4000 -C cca-mnn -z -c 40 -m 250 -v -y 200,4000 -u clustree,sil




##### Archived runs #########

# QC only

#Rscript /home/sonas/beegfs/my_scripts/pipeline/my_scAPA/seurat_for_scpolya/my_scRNA_pipeline_for_scpolya.R -i /home/sonas/beegfs/APA/scAPA/ureter10/peak_mat/l300iA10mm3/counts_mat_dir/ -o /home/sonas/beegfs/APA/scAPA/ureter10/peak_mat/l300iA10mm3/seurat_results/2022_09_26_u10_peak_mat_P0_200_Pn2000_qc_only/ -f 2022_09_26_u10_peak_mat_P0_200_Pn2000_qc_only -Q -z -c 30 -m 100 -v -q

#Rscript /home/sonas/beegfs/my_scripts/pipeline/my_scAPA/seurat_for_scpolya/my_scRNA_pipeline_for_scpolya.R -i /home/sonas/beegfs/APA/scAPA/ureter10/peak_mat/l300iA10mm3/counts_mat_dir/ -o /home/sonas/beegfs/APA/scAPA/ureter10/peak_mat/l300iA10mm3/seurat_results/2022_09_26_u10_peak_mat_P0_200_Pn2000/ -f 2022_09_26_u10_peak_mat_P0_200_Pn2000 -C cca-mnn -z -c 30 -m 100 -v -y 200,2000 


#Rscript /home/sonas/beegfs/my_scripts/pipeline/my_scAPA/seurat_for_scpolya/my_scRNA_pipeline_for_scpolya.R -i /home/sonas/beegfs/APA/scAPA/ureter10/peak_mat/l300iA10mm3/counts_mat_dir/ -o /home/sonas/beegfs/APA/scAPA/ureter10/peak_mat/l300iA10mm3/seurat_results/2022_09_26_u10_peak_mat_hvf3000/ -f 2022_09_26_u10_peak_mat_hvf3000 -C cca-mnn -z -c 30 -m 100 -v -y 3000 


#Rscript /home/sonas/beegfs/my_scripts/pipeline/my_scAPA/seurat_for_scpolya/my_scRNA_pipeline_for_scpolya.R -i /home/sonas/beegfs/APA/scAPA/ureter10/peak_mat/l300iA10mm3/counts_mat_dir/ -o /home/sonas/beegfs/APA/scAPA/ureter10/peak_mat/l300iA10mm3/seurat_results/2022_09_26_u10_peak_mat_P0_500_Pn4000/ -f 2022_09_26_u10_peak_mat_P0_500_Pn4000 -C cca-mnn -z -c 30 -m 100 -v -y 500,4000 


#Rscript /home/sonas/beegfs/my_scripts/pipeline/my_scAPA/seurat_for_scpolya/my_scRNA_pipeline_for_scpolya.R -i /home/sonas/beegfs/APA/scAPA/ureter10/peak_mat/l300iA10mm3/counts_mat_dir/ -o /home/sonas/beegfs/APA/scAPA/ureter10/peak_mat/l300iA10mm3/seurat_results/2022_09_26_u10_peak_mat_hvf6000/ -f 2022_09_26_u10_peak_mat_hvf6000 -C cca-mnn -z -c 30 -m 100 -v -y 6000 

Rscript /home/sonas/beegfs/my_scripts/pipeline/my_scAPA/seurat_for_scpolya/my_scRNA_pipeline_for_scpolya.R -i /home/sonas/beegfs/APA/scAPA/ureter10/filtered_gene_mat//home/sonas/beegfs/APA/scAPA/ureter10/filtered_gene_mat/counts_mat_dir/ -o /home/sonas/beegfs/APA/scAPA/ureter10/filtered_gene_mat/seurat_results/2022_09_28_u10_cell_mat_hvg2000/ -f 2022_09_28_u10_cell_mat_hvg2000 -C cca-mnn -z -c 30 -m 100 -v -y 2000






