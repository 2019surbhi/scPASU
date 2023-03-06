#!/usr/bin/bash

#SBATCH --mail-type=END
#SBATCH --mail-user=sonas@ccf.org
#SBATCH --job-name=saf
#SBATCH -n 20
#SBATCH --mem-per-cpu=3900

#SBATCH -o saf%J.out
#SBATCH -e saf_%J.err


### Load modules ###

module load R/4.1.0

ref_file_path=args[1]
out=args[2]
prefix=args[3]

ref_file='/home/sonas/beegfs/APA/scPASU/output/2_PeakRef/2e_TSS_filtered/u10_uro_l300iA10mm3_e100_updated_tss_filt.txt'
out='/home/sonas/beegfs/APA/scPASU/output/3_PeakCounts/'
prefix='u10_uro_l300iA10mm3_e100'

Rscript ${script_dir}2d_create.saf.R ${ref_file} ${out} ${prefix}
