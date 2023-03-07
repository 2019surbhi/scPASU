#!/usr/bin/bash

#SBATCH --mail-user=sonas@ccf.org
#SBATCH --job-name=bw_deepTools
#SBATCH -n 20 -p xtreme
#SBATCH --mem-per-cpu=3900

#SBATCH -o bam_to_bw_%J.out
#SBATCH -e bam_to_bw_%J.err

export LD_LIBRARY_PATH=/cm/local/apps/gcc/8.2.0/lib64/
module load gcc/8.2.0
module load samtools/1.9


#bam_dir='/home/sonas/beegfs/APA/scPASU/output/polyAFilter_optimization/l300A10mm0/'
#bam_dir='/home/sonas/beegfs/APA/scPASU/output/polyAFilter_optimization/l300iA6mm0/'
#bam_dir='/home/sonas/beegfs/APA/scPASU/output/1_process_bam/1c_filtered_bam/'

#bam_prefix='u10_uro_l300A10mm0'
#bam_prefix='u10_uro_l300iA6mm0'
#bam_prefix='u10_uro_merged'

out_dir='/home/sonas/beegfs/APA/scPASU/output/7_UCSC_uploads/bam_bw/bin10/'


norm=BPM
bin=10

bamCoverage --bam ${bam_dir}${bam_prefix}.bam -o ${out_dir}${bam_prefix}_${norm}_${bin}_minus.bw --samFlagInclude 16 --normalizeUsing ${norm} --binSize ${bin} -p 20

bamCoverage --bam ${bam_dir}${bam_prefix}.bam -o ${out_dir}${bam_prefix}_${norm}_${bin}_plus.bw --samFlagExclude 16 --normalizeUsing ${norm} --binSize ${bin} -p 20


### Run as loop ###

bam_dir='/home/sonas/beegfs/APA/scPASU/output/polyAFilter_optimization/U3/'

cd ${bam_dir}
bam=*.bam
prefix=${bam%.bam}

for bam_prefix in ${prefix[*]};do

bamCoverage --bam ${bam_dir}${bam_prefix}.bam -o ${out_dir}${bam_prefix}_${norm}_${bin}_minus.bw --samFlagInclude 16 --normalizeUsing ${norm} --binSize ${bin} -p 20

bamCoverage --bam ${bam_dir}${bam_prefix}.bam -o ${out_dir}${bam_prefix}_${norm}_${bin}_plus.bw --samFlagExclude 16 --normalizeUsing ${norm} --binSize ${bin} -p 20

done
