#!/usr/bin/bash

#SBATCH --mail-type=END
#SBATCH --mail-user=sonas@ccf.org
#SBATCH --job-name=10polya_counts_mat
#SBATCH -n 15
#SBATCH --array=1-10
#SBATCH --mem-per-cpu=3900

#SBATCH -o counts%J.out
#SBATCH -e counts%J.err

module load R/4.1.0 
module load python/3.9.5
module load bedtools
module load subset-bam/1.1.0

script_dir='/home/sonas/beegfs/APA/scPASU/scripts/'

### This requires .saf format ref which was created by the final peak ref creating script ###


### Initialize variables ###

# BAM #

#bam_dir='/home/sonas/beegfs/APA/scAPA/ureter10/dedup_ureter10_uro_filtered_BAM/'
#bam_dir='/home/sonas/beegfs/APA/scAPA/ureter10/dedup_ureter10_uro_filtered_BAM_l300_i_A10mm3/BAM/'
bam_dir='/home/sonas/beegfs/APA/scPASU/output/1_process_bam/1c_filtered_bam/'

sample_dir='/home/sonas/beegfs/APA/scPASU/input/'
export sample=$(cat ${sample_dir}ureter_samples.txt | nl -w1 -s ' ' | grep "^$SLURM_ARRAY_TASK_ID " | cut -f2 -d ' ')

#bam=${bam_dir}dedup_${sample}_uro_filtered.bam
#bam=${bam_dir}dedup_${sample}_uro_filtered_l300i_A10mm3.bam
bam=${bam_dir}dedup_${sample}_uro_filtered_l300iA10mm3.bam

# PEAK REF #

#peak_ref_dir='/home/sonas/beegfs/APA/scAPA/ureter10/peak_mat/final_ref/'
#peak_ref_dir='/home/sonas/beegfs/APA/scAPA/ureter10/rerun/peak_universe_w_extn/final_ref/'
peak_ref_dir='/home/sonas/beegfs/APA/scPASU/output/3_PeakCounts/'

#peak_ref=${peak_ref_dir}u10_uro_BAMfiltered_peak_universe.saf
#peak_ref=${peak_ref_dir}u10_uro_final_peak_ref_3polya_updated.saf
#peak_ref=${peak_ref_dir}u10_uro_final_peak_ref_5polya_updated.saf
#peak_ref=${peak_ref_dir}u10_uro_final_peak_ref_10polya_updated.saf
#peak_ref=${peak_ref_dir}u10_uro_final_peak_ref_l300i_A10mm3_3polya_updated.saf
#peak_ref=${peak_ref_dir}u10_uro_w_extn100_updated.saf
#peak_ref=${peak_ref_dir}u10_uro_w_extn100_updated_tss_filt.saf
peak_ref=${peak_ref_dir}u10_uro_l300iA10mm3_e100_updated.saf

# PEAK MAT DIR #

#peak_mat_dir='/home/sonas/beegfs/APA/scAPA/ureter10/peak_mat/l300iA10mm3/'
#peak_mat_dir1='/home/sonas/beegfs/APA/scAPA/ureter10/rerun/peak_mat/l300iA10mm3_extn100/'
#peak_mat_dir='/home/sonas/beegfs/APA/scAPA/ureter10/rerun/peak_mat/l300iA10mm3_extn100_tss_filt/'
peak_mat_dir='/home/sonas/beegfs/APA/scPASU/output/3_PeakCounts/'

### 2. Split barcode per cell ###

#mkdir ${peak_mat_dir}split_barcodes/${sample}
#Rscript ${script_dir}3a_split_bc.R ${peak_mat_dir} ${sample}


### 3. split bam per cell ###

#mkdir ${peak_mat_dir}split_bam/${sample}

subset_bam_dir=${peak_mat_dir}3b_SplitBAM/${sample}/
#bc_dir=${peak_mat_dir}3a_SplitBarcodes/${sample}/

#cd ${bc_dir}
#bc=*.tsv

#for barcode in ${bc[*]};do
#bc_name=${barcode%.tsv}
#subset-bam -b ${bam} -c ${bc_dir}${barcode} -o ${subset_bam_dir}${sample}_uro_${bc_name}.bam --cores 10 --bam-tag CB:Z
#done


### 4. Feature Counts ###

echo Generate peak count per cell
counts_dir=${peak_mat_dir}counts_dir/${sample}/
mkdir -p ${counts_dir}

cd ${subset_bam_dir}
b=*.bam

for bamfile in ${b[*]};do
pre=${bamfile%.bam}
Rscript ${script_dir}3b_feature_counts.R ${bamfile} ${peak_ref} ${counts_dir} ${pre} 15 no
done


### 5. Merge counts ###
echo Merge peak counts per sample

counts_mat_dir=${peak_mat_dir}counts_mat_dir/
mkdir -p ${counts_mat_dir}

prefix=${sample}_uro

Rscript ${script_dir}3c_merge_counts.R ${counts_dir}  ${prefix} ${counts_mat_dir}

