#!/usr/bin/bash

#SBATCH --mail-type=END
#SBATCH --mail-user=sonas@ccf.org
#SBATCH --job-name=peak_cov
#SBATCH -n 30
#SBATCH --mem-per-cpu=3900

#SBATCH -o peak_cov_%J.out
#SBATCH -e peak_cov_%J.err

module load R/4.1.0

### Initialize variables ###

#strand='minus'
#is_minus='TRUE'

strand='plus'
is_minus='FALSE'

cores=30
bam_dir='/home/sonas/beegfs/APA/scPASU/output/1_process_bam/1e_merged_bam/'
bam=${bam_dir}dedup_u10_uro_clean_filtered2_2_merged.bam

ref_dir='/home/sonas/beegfs/APA/scPASU/output/3_RefinePeakRef/3a_assign_TU/'
ref_saf=${ref_dir}u10_uro_filtered_100_peak_universe_updated.saf

outdir='/home/sonas/beegfs/APA/scPASU/output/3_RefinePeakRef/3b_PeakCoverage/'
fprefix=u10_uro_filtered_100

script_dir='/home/sonas/beegfs/APA/scPASU/scripts/'

# 1. Calculate peak counts

echo peak counts
Rscript ${script_dir}feature_counts.R ${bam} ${ref_saf} ${outdir} ${fprefix} ${cores} no

# 2. Calculate per TU peak coverage

counts_file=${outdir}${fprefix}_peak_count.rds
ref=${ref_dir}u10_uro_filtered_100_peak_universe_updated.txt

echo get peak coverage
Rscript ${script_dir}3b_peaks_coverage.R ${ref} ${counts_file} ${outdir}


# 3.  Plot peak coverage percent

echo peak coverage plot 
pcov_file=${outdir}u10_uro_filtered_100_peak_count_updated.rds
fprefix=u10_uro_filtered_100
w=11
h=8
x=100
y=5000

Rscript ${script_dir}3b_plot_peak_cov.R ${pcov_file} ${outdir} ${fprefix} ${w} ${h} ${x} ${y}

# 4. Filter peak count mat
echo filter peak ref

min_cov=10
out='/home/sonas/beegfs/APA/scPASU/output/3_RefinePeakRef/3c_updated_peaks/'
fprefix=u10_uro_filtered_100

Rscript ${script_dir}3b_filter_peak_ref.R ${ref} ${pcov_file} ${out} ${fprefix} ${min_cov} ${is_minus}

