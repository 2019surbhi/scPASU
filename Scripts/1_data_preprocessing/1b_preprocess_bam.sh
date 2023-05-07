#!/usr/bin/bash

#SBATCH --mail-type=END
#SBATCH --mail-user=sonas@ccf.org
#SBATCH --job-name=preprocess_bam
#SBATCH --array=1-10
#SBATCH -n 20
#SBATCH --mem-per-cpu=3900

#SBATCH -o preprocess_bam_%J.out
#SBATCH -e preprocess_bam_%J.err

export LD_LIBRARY_PATH=/cm/local/apps/gcc/8.2.0/lib64/
module load gcc/8.2.0
module load samtools/1.9
module load subset-bam/1.1.0
module load python/3.8.6 #To use umi_tools

# Input dir #
input_dir='/home/sonas/beegfs/APA/scPASU/input/'
barcodes=${input_dir}'/barcode/'
bam_dir=${input_dir}'/bam/'

# Output directories #
outdir='/home/sonas/beegfs/APA/scPASU/output/1_process_bam/'
subset_dir=${outdir}'/1a_subset_bam/'
dedup_dir=${outdir}'/1b_dedup_bam/'
clean_bam_dir=${outdir}'/1c_clean_bam/'
filter_dir=${outdir}'/1d_filtered_bam/'

export samples_file=${input_dir}'ureter_samples.txt'
sample=$(cat ${samples_file} | nl -w1 -s ' ' | grep "^$SLURM_ARRAY_TASK_ID " | cut -f2 -d ' ')

## 1. Subset ##

subset-bam -b ${bam_dir}${sample}.bam -c ${barcodes}${sample}_uro_barcodes.tsv -o ${subset_dir}${sample}_uro.bam --log-level debug
samtools index ${subset_dir}${sample}_uro.bam

## 2. Dedup ##

umi_tools dedup -I ${subset_dir}${sample}_uro.bam -S ${dedup_dir}dedup_${sample}_uro.bam --method=unique --extract-umi-method=tag --umi-tag=UB --cell-tag=CB
samtools index ${dedup_dir}dedup_${sample}_uro.bam

## 3. Cleanup BAM - only uniquely mapped reads retained ##

#samtools view -b -q 1 ${dedup_dir}dedup_${sample}_uro.bam > ${clean_bam_dir}dedup_uniq_${sample}_uro.bam
#samtools index ${clean_bam_dir}dedup_uniq_${sample}_uro.bam

## 4. Filter genomicA reads ##

WORKING_DIR='/home/sonas/beegfs/tools/polyAfilter/'

FASTA="/home/sonas/beegfs/ref/refdata-gex-GRCh38-2020-A/fasta/genome.fa"
GTF="/home/sonas/beegfs/ref/refdata-gex-GRCh38-2020-A/genes/genes.gtf"

BAM=${clean_bam_dir}dedup_uniq_${sample}_uro.bam

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

