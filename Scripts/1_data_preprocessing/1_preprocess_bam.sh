#!/usr/bin/bash

#SBATCH --mail-type=END
#SBATCH --mail-user=sonas@ccf.org
#SBATCH --job-name=preprocess_bam
#SBATCH --array=1-10
#SBATCH -n 20
#SBATCH --mem-per-cpu=3900

#SBATCH -o preprocess_bam_ureter%J.out
#SBATCH -e preprocess_bam_ureter%J.err

export LD_LIBRARY_PATH=/cm/local/apps/gcc/8.2.0/lib64/
module load gcc/8.2.0
module load samtools/1.9
module load subset-bam/1.1.0
module load python/3.8.6 #To use umi_tools

barcodes='/home/sonas/beegfs/APA/scAPA/ureter10/uro_barcodes/'
bam_dir='/home/sonas/beegfs/APA/scAPA/ureter10/ureter_bam/'
subset_dir='/home/sonas/beegfs/APA/scAPA/ureter10/u10_uro_bam/'
dedup_dir='/home/sonas/beegfs/APA/scAPA/ureter10/dedup_u10_uro_bam/'
filter_dir='/home/sonas/beegfs/APA/scAPA/ureter10/dedup_u10_uro_filterd_bam/'

export samples_file='/home/sonas/beegfs/APA/scAPA/ureter10/ureter.txt'



sample=$(cat ${samples_file} | nl -w1 -s ' ' | grep "^$SLURM_ARRAY_TASK_ID " | cut -f2 -d ' ')

## 1. Subset ##

subset-bam -b ${bam_dir}${sample}.bam -c ${barcodes}${sample}_uro_barcodes.tsv -o ${subset_dir}${sample}_uro.bam --log-level debug

samtools index ${subset_dir}${sample}_uro.bam

## 2. Dedup ##

umi_tools dedup -I ${subset_dir}${sample}_uro.bam -S ${dedup_dir}dedup_${sample}_uro.bam --method=unique --extract-umi-method=tag --umi-tag=UB --cell-tag=CB

samtools index ${dedup_dir}dedup_${sample}_uro.bam

## 3. Filter genomicA reads ##

WORKING_DIR='/home/sonas/beegfs/tools/polyAfilter/'

FASTA="/home/sonas/beegfs/ref/refdata-gex-GRCh38-2020-A/fasta/genome.fa"
GTF="/home/sonas/beegfs/ref/refdata-gex-GRCh38-2020-A/genes/genes.gtf"

BAM=${dedup_dir}dedup_${sample}_uro.bam

DB=${filter_dir}"gtfdb.db"
TRANS=${filter_dir}${sample}"_trans.pkl" # Note this is a cache file so define with sample name if running in parallel

NTHREADS=20

#SCRIPT
cd $WORKING_DIR

# Create the GTF database
python polyAfilter.py createDB -v $GTF $DB

# Create the TRANS file
python polyAfilter.py createTRANS $DB $BAM $TRANS

# Filter the BAM file
python polyAfilter.py BAMfilter -m 1 -v -p $NTHREADS 300 10 $BAM $FASTA $TRANS -o ${filter_dir}dedup_${sample}_uro_filtered.bam

samtools index ${filter_dir}dedup_${sample}_uro_filtered.bam

