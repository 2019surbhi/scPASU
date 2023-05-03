#!/usr/bin/bash

#SBATCH --mail-type=END
#SBATCH --mail-user=sonas@ccf.org
#SBATCH --job-name=split_bam
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=3900

#SBATCH -o split_bam_m.out
#SBATCH -e split_bam_m.err

export LD_LIBRARY_PATH=/cm/local/apps/gcc/8.2.0/lib64/
module load gcc/8.2.0
module load samtools/1.9
module load subset-bam/1.1.0

bam_dir='/home/sonas/beegfs/APA/scPASU/output/1_process_bam/1d_merged_bam/'
#out_dir='/home/sonas/beegfs/APA/scPASU/output/1_process_bam/1d_merged_bam/new/'

#out_dir='/home/sonas/beegfs/APA/scPASU/output/1_process_bam/1d_merged_bam/new2/'

prefix='u10_uro'
bam=${bam_dir}${prefix}_merged.bam

### Single flag ###

#cd  ${out_dir}

#samtools view -b -F 0x10 ${bam} > ${prefix}_plus.bam
#samtools index ${prefix}_plus.bam

#samtools view -b -f 0x10 ${bam} > ${prefix}_minus.bam
#samtools index ${prefix}_minus.bam


#out_dir='/home/sonas/beegfs/APA/scPASU/output/1_process_bam/1d_merged_bam/new3/'

#cd ${out_dir}

#samtools view -b -h -F 20 ${bam} > ${prefix}_plus.bam
#samtools index ${prefix}_plus.bam

#samtools view -b -h -f 16 ${bam} > ${prefix}_minus.bam
#samtools index ${prefix}_minus.bam

### Cleanup + split ###

out_dir='/home/sonas/beegfs/APA/scPASU/output/1_process_bam/1d_merged_bam/new4/'

cd ${out_dir}

samtools view -b -F 4 ${bam} > mapped.bam
samtools index mapped.bam

samtools view -b -q 1 mapped.bam > uniq_mapped.bam
samtools index uniq_mapped.bam


new_bam=unique.mapped.bam

#samtools view -b -h -F 20 ${new_bam} > ${prefix}_plus.bam
#samtools index ${prefix}_plus.bam

#samtools view -b -h -f 16 ${new_bam} > ${prefix}_minus.bam
#samtools index ${prefix}_minus.bam



## Split into plus strand ##

cd ${out_dir}

#samtools view -b -f 128 -F 16 ${bam} > fwd1.bam
#samtools index fwd1.bam

#samtools view -b -f 80 ${bam} > fwd2.bam

#samtools view -b -F 16 ${bam} > fwd2.bam
#samtools index fwd2.bam

# Combine alignments that originate on the forward strand.

#samtools merge -f ${prefix}_plus.bam fwd1.bam fwd2.bam
#samtools index ${prefix}_plus.bam

## Split into minus strand ##

#cd ${out_dir}
#samtools view -b -f 144 ${bam} > rev1.bam
#samtools index rev1.bam

#samtools view -b -f 64 -F 16 ${bam} > rev2.bam

#samtools view -b -f 16 ${bam} > rev2.bam
#samtools index rev2.bam

# Combine alignments that originate on the reverse strand.

#samtools merge -f ${prefix}_minus.bam rev1.bam rev2.bam
#samtools index ${prefix}_minus.bam

#### Try again ####

# Ref: https://www.seqanswers.com/forum/bioinformatics/bioinformatics-aa/24884-strand-specific-bam-files#post208995


# Plus #

#cd ${out_dir}

#samtools view -b -f 99 ${bam} > fwd1.bam
#samtools index fwd1.bam

#samtools view -b -f 147 ${bam} > fwd2.bam
#samtools index fwd2.bam

# Combine alignments that originate on the forward strand.

#samtools merge -f ${prefix}_plus.bam fwd1.bam fwd2.bam
#samtools index ${prefix}_plus.bam
# Minus #

#cd ${out_dir}
#samtools view -b -f 83 ${bam} > rev1.bam
#samtools index rev1.bam

#samtools view -b -f 163 ${bam} > rev2.bam
#samtools index rev2.bam

# Combine alignments that originate on the reverse strand.

#samtools merge -f ${prefix}_minus.bam rev1.bam rev2.bam
#samtools index ${prefix}_minus.bam
