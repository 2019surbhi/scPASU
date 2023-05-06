#!/usr/bin/bash

#SBATCH --mail-user=sonas@ccf.org
#SBATCH --job-name=seurat_for_scpolya
#SBATCH -n 40
#SBATCH --mem-per-cpu=3900

#SBATCH -o seurat_for_scpolya.out
#SBATCH -e seurat_for_scpolya.err

module load R

apa_file='/home/sonas/beegfs/APA/scAPA/ureter10/misc/apa_genes.txt'
Rscript my_scRNA_pipeline_for_scpolya.R -i /home/sonas/beegfs/APA/scPASU/output/3_PeakCounts/3c_CreatePeakMat/counts_mat_dir/ -o ${peak_mat_clus_dir}2023_03_06_P0_200_Pn_4000/ -f 2023_03_06_P0_200_Pn_4000 -C cca-mnn -u clustree,sil -z -c 30 -m 100 -v -d 1:50 -y 200,4000

