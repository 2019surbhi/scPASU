#!/usr/bin/bash

#SBATCH --mail-type=END
#SBATCH --mail-user=sonas@ccf.org
#SBATCH --job-name=motif
#SBATCH -n 35 -p xtreme
#SBATCH --mem-per-cpu=7900

#SBATCH -o motif_pr_plus_us.out
#SBATCH -e motif_pr_plus_us.err


module load meme/5.4.1

motif_dir='/home/sonas/beegfs/APA/scPASU/output/2.5_PeakRefAssessment/motif_analysis/'
fa_file=${motif_dir}pr_us_fasta.fa
fa_file2=${motif_dir}pr_fasta.fa
rna_alpha=${motif_dir}rna_alpa.fa


# Run it for upstream region
#streme --p ${fa_file} --o ${motif_dir}pr_upstream/ --w 6

# Run it for PR+upstream region
streme --p ${fa_file2} --o ${motif_dir}pr_plus_upstream/ --w 6




