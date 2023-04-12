#!/usr/bin/bash

#SBATCH --mail-type=END
#SBATCH --mail-user=sonas@ccf.org
#SBATCH --job-name=create_peak_ref_BAM
#SBATCH -n 30 -p xtreme
#SBATCH --mem-per-cpu=3900

#SBATCH -o peak_ref_filtered_BAM.out
#SBATCH -e peak_ref_filtered_BAM.err


### Load modules ###

# conda activate scPASU

peak_ref='/home/sonas/beegfs/APA/scPASU/output/2_PeakRef/2e_TSS_filtered/u10_uro_l300iA10mm3_e100_updated_tss_filt_final.bed'
BAM='/home/sonas/beegfs/APA/scPASU/output/1_process_bam/1d_merged_bam/u10_uro_merged.bam'
out='/home/sonas/beegfs/APA/scPASU/output/7_UCSC_uploads/peak_filtered_bam/'

bedtools intersect -abam ${BAM} -b ${peak_ref} -wa > ${out}u10_uro_peak_filtered.bam

bam_dir='/home/sonas/beegfs/APA/scPASU/output/3_PeakCounts/3b_SplitBAM_per_clus/'
bam_pre='2021_04_01_ureter10_uro_PC50_res0.2_cluster'
peak_ref='/home/sonas/beegfs/APA/scPASU/output/2_PeakRef/2e_TSS_filtered/u10_uro_l300iA10mm3_e100_updated_tss_filt_final.bed'
out='/home/sonas/beegfs/APA/scPASU/output/7_UCSC_uploads/peak_filtered_bam/'

clus=$SLURM_ARRAY_TASK_ID

bedtools intersect -abam ${bam_dir}${bam_pre}${clus}.bam -b ${peak_ref} -wa > ${out}cluster${clus}_peak_filtered.bam
