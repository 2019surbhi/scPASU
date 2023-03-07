#!/usr/bin/bash

#SBATCH --mail-user=sonas@ccf.org
#SBATCH --job-name=seurat_for_scpolya
#SBATCH -n 40
#SBATCH --mem-per-cpu=3900

#SBATCH -o seurat_for_scpolya.out
#SBATCH -e seurat_for_scpolya.err

module load R

apa_file='/home/sonas/beegfs/APA/scAPA/ureter10/misc/apa_genes.txt'
Rscript /home/sonas/beegfs/APA/scAPA/ureter10/misc/my_scRNA_pipeline_for_scpolya.R -i /home/sonas/beegfs/APA/scAPA/ureter10/rerun/peak_mat/l300iA10mm3_extn100_tss_filt/counts_mat_dir/ -o /home/sonas/beegfs/APA/scAPA/ureter10/rerun/peak_mat/l300iA10mm3_extn100_tss_filt/seurat_results/2022_12_25_u10_peak_mat_extn100_tss_filt_P0_200_Pn_APA/ -f 2022_12_25_u10_peak_mat_extn100_tss_filt_P0_200_Pn_APA -C cca-mnn -z -c 40 -m 250 -v -y 200,2000 -Y ${apa_file}


