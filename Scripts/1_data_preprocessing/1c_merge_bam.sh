#!/usr/bin/bash

#SBATCH --mail-type=END
#SBATCH --mail-user=sonas@ccf.org
#SBATCH --job-name=hg38_prepare_bam
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=3900

#SBATCH -o prepare_bam_ureterhg38.out
#SBATCH -e prepare_bam_ureterhg38.err

export LD_LIBRARY_PATH=/cm/local/apps/gcc/8.2.0/lib64/
module load gcc/8.2.0
module load samtools/1.9

# Merge unfiltered BAM (BEFORE genomicA filtering)
bam_dir='/home/sonas/beegfs/APA/scPASU/output/1_process_bam/1c_clean_bam/'
prefix='dedup_u10_uro_clean'

# Merge unfiltered BAM (AFTER to genomicA filtering)
bam_dir='/home/sonas/beegfs/APA/scPASU/output/1_process_bam/1d_filtered_bam/'
prefix='dedup_u10_uro_clean_filter'


out_dir='/home/sonas/beegfs/APA/scPASU/output/1_process_bam/1e_merged_bam/'

### Merge filtered BAMs ###

cd ${bam_dir}
bam=$(ls --hide=*.bai)

## a. Merge  ##

cd ${bam_dir}

samtools merge ${out_dir}${prefix}_merged.bam ${bam}
samtools index ${out_dir}${prefix}_merged.bam

## b. Split by strand ##

samtools view -F 16 -b ${out_dir}${prefix}_merged.bam > ${out_dir}${prefix}_plus.bam
samtools view -f 16 -b ${out_dir}${prefix}_merged.bam > ${out_dir}${prefix}_minus.bam

samtools index ${out_dir}${prefix}_plus.bam
samtools index ${out_dir}${prefix}_minus.bam
