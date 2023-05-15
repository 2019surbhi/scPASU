#!/usr/bin/bash

#SBATCH --mail-type=END
#SBATCH --mail-user=sonas@ccf.org
#SBATCH --job-name=update_peak_ref
#SBATCH -n 30
#SBATCH --mem-per-cpu=3900

#SBATCH -o update_peak_ref_%J.out
#SBATCH -e update_peak_ref_%J.err

module load R/4.1.0

### Initialize variables ###

cores=30
script_dir='/home/sonas/beegfs/APA/scPASU/scripts/'
fprefix='u10_uro'

ref_dir='/home/sonas/beegfs/APA/scPASU/output/3_RefinePeakRef/3a_assign_TU/'
ref=${ref_dir}${fprefix}_peak_universe_updated.txt

outdir='/home/sonas/beegfs/APA/scPASU/output/3_RefinePeakRef/'

mkdir -p ${outdir}

# 1a. Filter by PR length

echo filter by PR width

maxwidth=1000
out1=${outdir}3b_Filter_by_PR_width/
mkdir -p ${out1}

Rscript ${script_dir}3b_filter_by_pr_width.R ${ref} ${out1} ${fprefix} ${maxwidth}

# 1b. Assign TU - to update peak ref

echo  2nd TU assignment

dir='/home/sonas/beegfs/APA/scPASU/'
peak_ref_dir='none'
ref2=${out1}${fprefix}_filtered_peak_ref.txt
fprefix2=u10_uro_filtered_tu_assigned

Rscript ${script_dir}3a_assign_tu.R 30 ${dir} ${peak_ref_dir} ${fprefix2} ${out1} ${gtf_file} ${ref2}

# 2a. Calculate peak counts

echo get peak counts

bam_dir='/home/sonas/beegfs/APA/scPASU/output/1_process_bam/1e_merged_bam/'
bam=${bam_dir}dedup_u10_uro_clean_filtered2_2_merged.bam

ref_saf=${out1}${fprefix2}_peak_universe_updated.saf
out2=${outdir}'3c_PeakCoverage/'

Rscript ${script_dir}feature_counts.R ${bam} ${ref_saf} ${out2} ${fprefix} ${cores} no

# 2b. Calculate per TU peak coverage

echo get peak coverage

counts_file=${out2}${fprefix}_peak_count.rds
ref3=${out1}${fprefix2}_peak_universe_updated.txt

Rscript ${script_dir}3b1_peaks_coverage.R ${ref3} ${counts_file} ${out2}

# 2c.  Plot peak coverage percent

echo peak coverage plot 
pcov_file=${out2}${fprefix}_peak_count_updated.rds
fprefix=u10_uro_filtered_100
w=11
h=8
x=100
y=5000

Rscript ${script_dir}3b2_plot_peak_cov.R ${pcov_file} ${out2} ${fprefix} ${w} ${h} ${x} ${y}

# 3. Filter peak count mat
echo filter peak ref

out3=${outdir}'3d_final_ref/'
mkdir -p ${out3}

min_cov=10

fprefix_final=${fprefix}_final
pcov_file=${out2}${fprefix}_peak_count_updated.rds

Rscript ${script_dir}3c_filter_by_peak_cov.R ${ref} ${pcov_file} ${out} ${fprefix} ${min_cov} ${is_minus}

