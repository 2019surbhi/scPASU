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





