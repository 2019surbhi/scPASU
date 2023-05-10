#!/usr/bin/bash

#SBATCH --mail-type=END
#SBATCH --mail-user=sonas@ccf.org
#SBATCH --job-name=final_peak_ref
#SBATCH -n 20 -p xtreme
#SBATCH --mem-per-cpu=3900

#SBATCH -o peak_cov%J.out
#SBATCH -e peak_cov%J.err


### Load modules ###

module load R/4.1.0

### Initialize variables ###

strand='minus'
strand='plus'

cores=30
bam_dir='/home/sonas/beegfs/APA/scPASU/output/1_process_bam/1e_merged_bam/'
bam=dedup_u10_uro_clean_filtered2_${strand}.bam

ref='/home/sonas/beegfs/APA/scPASU/output/3_RefinePeakRef/3a_assign_TU/u10_uro_peak_universe_updated.saf'
outdir='/home/sonas/beegfs/APA/scPASU/output/3_RefinePeakRef/3b_PeakCoverage'
fprefix='u10_uro'

script_dir='/home/sonas/beegfs/APA/scPASU/scripts/'
Rscript ${script_dir}feature_counts.R ${bam} ${ref} ${outdir} ${fprefix} ${cores} no



