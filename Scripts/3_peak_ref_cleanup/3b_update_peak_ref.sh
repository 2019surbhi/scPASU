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
is_minus='TRUE'

#strand='plus'
#is_minus='FALSE'

cores=30
bam_dir='/home/sonas/beegfs/APA/scPASU/output/1_process_bam/1e_merged_bam/'
bam=${bam_dir}dedup_u10_uro_clean_filtered_${strand}.bam

ref_dir='/home/sonas/beegfs/APA/scPASU/output/3_RefinePeakRef/3a_assign_TU/'
#ref=${ref_dir}u10_uro_peak_universe_updated.saf
ref=${ref_dir}u10_uro_peak_universe_${strand}_updated.saf

outdir='/home/sonas/beegfs/APA/scPASU/output/3_RefinePeakRef/3b_PeakCoverage2/'
fprefix=u10_uro_${strand}

script_dir='/home/sonas/beegfs/APA/scPASU/scripts/'

# 1. Calculate peak counts
#Rscript ${script_dir}feature_counts.R ${bam} ${ref} ${outdir} ${fprefix} ${cores} no

# 2. Calculate per TU peak coverage

counts_file=${outdir}${fprefix}_peak_count.rds

#Rscript ${script_dir}3b_peaks_coverage.R ${ref} ${counts_file} ${outdir}

# 3.  Plot peak coverage percent

pcov_file=${outdir}u10_uro_${strand}_peak_count_updated.rds
fprefix=u10_uro_${strand}
#w=11
#h=8
#x=100
#y=5000

#Rscript ${script_dir}3b_plot_peak_cov.R ${pcov_file} ${outdir} ${fprefix} ${w} ${h} ${x} ${y}

# 4. Filter peak count mat
min_cov=5

Rscript ${script_dir}3b_filter_peak_count.R ${pcov_file} ${outdir} ${fprefix} ${min_cov} ${is_minus}


