#!/usr/bin/bash

#SBATCH --mail-type=END
#SBATCH --mail-user=sonas@ccf.org
#SBATCH --job-name=cleaned_gene_mat
#SBATCH --array=1-10
#SBATCH -n 32
#SBATCH --mem-per-cpu=3900

#SBATCH -o cleaned_gene_mat%J.out
#SBATCH -e cleaned_gene_mat%J.err

module load cellranger/4.0.0
module load R/4.1.0 
module load STAR/2.7.8a


export sample=$(cat /home/sonas/beegfs/APA/scAPA/ureter10/ureter_samples.txt | nl -w1 -s ' ' | grep "^$SLURM_ARRAY_TASK_ID " | cut -f2 -d ' ')

#sample='U_199749'

## (A) Cellranger workflow ##

## Step 1: Covert BAM to fastq ##

#bam_dir='/home/sonas/beegfs/APA/scAPA/ureter10/dedup_ureter10_uro_filtered_BAM_l300_i_A10mm3/BAM/'
#bam=${bam_dir}dedup_${sample}_uro_filtered_l300i_A10mm3.bam

# Make sure the final folder/dir in the fastq_dir path doesn't exist else the script runs into "already exists" error!

#fastq_dir='/home/sonas/beegfs/APA/scAPA/ureter10/filtered_gene_mat/fastq/'

#cellranger bamtofastq --nthreads=20  ${bam} ${fastq_dir}/

## Step 2: Rename fastq files (if all have bamtofastq prefix) ##

#cd ${fastq_dir}
#s=*
#for sample in ${s[*]};do
#cd ${sample}
#files=*
#for file in ${files[*]};do
#new=${file/bamtofastq/$sample}
#mv ${file} ${new}
#done
#cd ..
#done


## Step 3: Map using cellranger ##

#memory=120
#cores=32
#ref_dir='/home/sonas/beegfs/ref/refdata-gex-GRCh38-2020-A/'
#counts_dir='/home/sonas/beegfs/APA/scAPA/ureter10/filtered_gene_mat/filtered_counts/'

#cd ${counts_dir}

## cellranger was over-filtering cells and drastically filtering the genes, so I added additional parameters ##

#echo "cellranger count --id ${sample} --fastqs=${fastq_dir}${sample} --sample=${sample} --localcores=${cores} --localmem=${memory} --transcriptome=${ref_dir} --no-target-umi-filter --target-panel=/cm/shared/apps/cellranger/4.0.0/target_panels/gene_signature_v1.0_GRCh38-2020-A.target_panel.csv"


#cellranger count --id ${sample} --fastqs=${fastq_dir}${sample} --sample=${sample} --localcores=${cores} --localmem=${memory} --transcriptome=${ref_dir} --no-target-umi-filter --target-panel=/cm/shared/apps/cellranger/4.0.0/target_panels/gene_signature_v1.0_GRCh38-2020-A.target_panel.csv



## (B) STAR solo workflow ##

#bam_dir='/home/sonas/beegfs/APA/scAPA/ureter10/dedup_ureter10_uro_filtered_BAM_l300_i_A10mm3/BAM/'

#bam=${bam_dir}dedup_${sample}_uro_filtered_l300i_A10mm3.bam

#STAR --runThreadN 30 \
#--runMode alignReads \
#--soloType CB_UMI_Simple \
#--genomeDir /home/sonas/beegfs/ref/STAR/refdata-gex-GRCh38-2020-A_star/ \
#--genomeFastaFiles /home/sonas/beegfs/ref/refdata-gex-GRCh38-2020-A/fasta/genome.fa \
#--sjdbGTFfile /home/sonas/beegfs/ref/refdata-gex-GRCh38-2020-A/genes/genes.gtf \
#--readFilesType SAM SE \
#--readFilesCommand samtools view -F 0x100 \
#--soloInputSAMattrBarcodeSeq CR UR \
#--soloInputSAMattrBarcodeQual CY UY \
#--readFilesSAMattrKeep None \
#--soloFeatures GeneFull \
#--soloCBwhitelist /home/sonas/beegfs/ref/STAR/3M-february-2018.txt \
#--outFileNamePrefix /home/sonas/beegfs/APA/scAPA/ureter10/filtered_gene_mat/starSolo_filtered_counts/${sample}/ \
#--readFilesIn ${bam} \


## (C) Use R subread to just get a count matrix ##

script_dir='/home/sonas/beegfs/my_scripts/pipeline/my_scAPA/peak_mat/'
peak_mat_dir='/home/sonas/beegfs/APA/scAPA/ureter10/peak_mat/l300iA10mm3/'
subset_bam_dir=${peak_mat_dir}split_bam/${sample}/


mkdir /home/sonas/beegfs/APA/scAPA/ureter10/filtered_gene_mat/counts_dir/${sample}
counts_dir=/home/sonas/beegfs/APA/scAPA/ureter10/filtered_gene_mat/counts_dir/${sample}/
gtf='/home/sonas/beegfs/ref/refdata-gex-GRCh38-2020-A/genes/genes.gtf'


cd ${subset_bam_dir}
b=*.bam

for bamfile in ${b[*]};do
pre=${bamfile%.bam}
Rscript ${script_dir}featurecounts.R ${bamfile} ${gtf} ${counts_dir} ${pre} 15 yes
done

### 5. Merge counts ###

counts_mat_dir='/home/sonas/beegfs/APA/scAPA/ureter10/filtered_gene_mat/counts_mat_dir/'

prefix=${sample}_uro

Rscript ${script_dir}merge_counts.R ${counts_dir}  ${prefix} ${counts_mat_dir}



