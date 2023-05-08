#!/usr/bin/bash

#SBATCH --mail-type=END
#SBATCH --mail-user=sonas@ccf.org
#SBATCH --job-name=plus_MACS2_peak_calling
#SBATCH --cpus-per-task=30
#SBATCH --mem-per-cpu=3900

#SBATCH -o macs2_u10_uro_p.out
#SBATCH -e macs2_u10_uro_p.err

#module load R/4.1.0

module load samtools


bam_dir='/home/sonas/beegfs/APA/scPASU/output/1_process_bam/1e_merged_bam/'
out_dir='/home/sonas/beegfs/APA/scPASU/output/2_PeakRef/2a_peaks/'

strand='plus'
strand='minus'

bam_prefix=dedup_u10_uro_clean_filtered_${strand}

# Run MACS2

cd ${bam_dir}

macs2 callpeak \
	-t ${bam_dir}${bam_prefix}.bam \
	--keep-dup all \
	--gsize hs \
	--format BAM \
	--nomodel \
        --extsize 200 \
        --shift -100 \
	--outdir ${out_dir} \
	--name ${bam_prefix} \
	-B \
	--verbose 2 
