#!/usr/bin/bash

#SBATCH --mail-type=END
#SBATCH --mail-user=sonas@ccf.org
#SBATCH --job-name=assign_tu
#SBATCH -n 30 -p xtreme
#SBATCH --mem-per-cpu=3900

#SBATCH -o assign_tu.out
#SBATCH -e assign_tu.err


### Load modules ###

module load R/4.1.0
module load python/3.9.5

module load bedtools/2.29.0

### Initialize variables ###


#extn=25
#extn=50
#extn=75
extn=100
#extn=125
#extn=149

# dir='/home/sonas/beegfs/APA/scPASU/'
#peaks_ref_dir='/home/sonas/beegfs/APA/scPASU/output/2_PeakRef/2c_intersect/e'${extn}'/'
#file_prefix='u10_uro'
#gtf_file='/home/sonas/beegfs/ref/refdata-gex-GRCh38-2020-A/genes/genes.gtf'
#out_dir=${dir}'output/3_RefinePeakRef/3a_assign_TU/'

#### Updated path for Ninh and Justin ###

#Please change the you_dir path to where you place the scPASU dir 
#Create `input' and 'output' dir within the scPASU dir. 
#Then Download the Github 'Scripts' dir withing scPASU folder [the paths in the R script should work then]
dir='/path/to/your_dir/scPASU/'

# Create 3_RefinePeakRef/3a_assign_TU/e100/ (Go to scPASU/output directory and use mkdir -p 3_RefinePeakRef/3a_assign_TU/
# Then download and place e100 within 3a_assign_TU
peaks_ref_dir=${dir}/output/2_PeakRef/2c_intersect/e100/

# Change the path to where you have gtf file in your reference folder
gtf_file='/home/sonas/beegfs/ref/refdata-gex-GRCh38-2020-A/genes/genes.gtf'

out_dir=${dir}/output/3_RefinePeakRef/3a_assign_TU/'

file_prefix='u10_uro'

Rscript ${dir}/scripts/2_assign_tu.R 30 ${dir} ${peaks_ref_dir} ${file_prefix} ${out_dir} ${gtf_file}

# You should get the following files from Box and add them within your ${dir}/scPASU/input/ folder:
# 1. GRCh38_gene_biotypes.csv
# 2. genes.rds (optional)
