#!/usr/bin/bash

#SBATCH --mail-type=END
#SBATCH --mail-user=sonas@ccf.org
#SBATCH --job-name=10polya_counts_mat
#SBATCH -n 30 -p xtreme
#SBATCH --mem-per-cpu=7900

#SBATCH -o counts%J.out
#SBATCH -e counts%J.err

module load R/4.1.0
module load python/3.9.5
module load bedtools
module load subset-bam/1.1.0

script_dir='/home/sonas/beegfs/APA/scPASU/scripts/'

sample_dir='/home/sonas/beegfs/APA/scPASU/input/'
peak_count_dir='/home/sonas/beegfs/APA/scPASU/output/3_PeakCounts/'
peak_ref=${peak_count_dir}u10_uro_l300iA10mm3_e100_updated.saf

bc_dir=${peak_count_dir}3a_SplitBarcodes_per_clus/
subset_bam_dir=${peak_count_dir}3b_SplitBAM_per_clus/


### 3b Subset BAM by cluster ###

#cd ${bc_dir}

#bc=*.tsv

#for barcode in ${bc[*]};do
#bc_name=${barcode%.tsv}
#subset-bam -b ${merged_bam} -c ${bc_dir}${barcode} -o ${subset_bam_dir}${bc_name}.bam --cores 30 --bam-tag CB:Z
#done

### 4b. Feature counts per clus ###

counts_dir=${peak_count_dir}3c_CreatePeakMat_per_clus/counts_dir/
mkdir -p ${counts_dir}

logdir=${peak_count_dir}3c_CreatePeakMat_per_clus/count_logs/
mkdir -p logdir

cd ${subset_bam_dir}
b=*.bam

for bamfile in ${b[*]};do
pre=${bamfile%.bam}
Rscript ${script_dir}3b_feature_counts.R ${bamfile} ${peak_ref} ${counts_dir} ${pre} 30 no ${logdir}
done


### 5b. Merge per clus counts ###
echo Merge peak counts per clus

counts_mat_dir=${peak_mat_dir}3c_CreatePeakMat_per_clus/counts_mat_dir/
mkdir -p ${counts_mat_dir}

prefix=${sample}_uro

Rscript ${script_dir}3c_merge_counts.R ${counts_dir}  ${prefix} ${counts_mat_dir}

