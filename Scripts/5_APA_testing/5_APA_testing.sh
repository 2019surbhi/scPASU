#!/usr/bin/bash

#SBATCH --mail-type=END
#SBATCH --mail-user=sonas@ccf.org
#SBATCH --job-name=APA_testing
#SBATCH -n 30 -p xtreme
#SBATCH --mem-per-cpu=3900

#SBATCH -o APA_testing.out
#SBATCH -e APA_testing.err


### Load modules ###

module load R/4.1.0

script_dir='/home/sonas/beegfs/APA/scPASU/scripts/'

### 1. Merge per sample counts to a single matrix ###

counts_dir='/home/sonas/beegfs/APA/scPASU/output/4_PeakCounts/4d_merged_PeakMat/'
fprefix='u101_uro'
out='/home/sonas/beegfs/APA/scPASU/output/5_APA_testing/inputs/'

Rscript ${script_dir}5a_merge_per_sample_counts.R ${counts_dir} ${fprefix} ${out}

### 2. APA testing ###

inputdir='/home/sonas/beegfs/APA/scPASU/output/5_APA_testing/inputs/'
outdir='/home/sonas/beegfs/APA/scPASU/output/5_APA_testing/outputs/'
counts_file='u10_uro_counts.txt'
meta_file='2021_04_01_ureter10_uro_PC50_res0.2_meta.xlsx'
peak_ref_file='/home/sonas/beegfs/APA/scPASU/output/5_APA_testing/inputs/u10_uro_APA_testing_peak_ref.txt'

mkdir -p outdir

Rscript ${script_dir}5_APA_testing_on_peak_mat.R ${inputdir} ${outdir} ${counts_file} ${meta_file} ${peak_ref_file}




