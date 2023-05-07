#!/usr/bin/bash

#SBATCH --mail-type=END
#SBATCH --mail-user=sonas@ccf.org
#SBATCH --job-name=polyA
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=3900

#SBATCH -o sam_to_polyA.out
#SBATCH -e sam_to_polyA.err

# conda activate my_scAPA

export LD_LIBRARY_PATH=/cm/local/apps/gcc/8.2.0/lib64/

module load gcc/8.2.0
module load samtools/1.9
#module load perl/5.34.0

PATH='/home/sonas/beegfs/miniconda3/envs/perl/lib/perl5/site_perl/5.22.0/Bio/DB/':$PATH

strand_dir='/home/sonas/beegfs/cellranger/mapped_data/ureter10/ureter10_uro/scAPA_1/temp/dedup_stranded_uro_bam/'
polya_dir='/home/sonas/beegfs/cellranger/mapped_data/ureter10/ureter10_uro/scAPA_1/temp/dedup_stranded_uro_polyA_sites/'
pola_genomicA_filtered_dir='/home/sonas/beegfs/cellranger/mapped_data/ureter10/ureter10_uro/scAPA_1/temp/pola_genomicA_filtered/'


cd ${strand_dir}
f=*minus.bam
for bam_file in ${f[*]};do

sample=${bam_file%_dedup_uro_minus.bam}

samtools view ${bam_file} |/home/sonas/beegfs/apa_tinglab/01_polyAseq/files_modified_for_scAPA/samToPolyA.pl --minClipped=10 --minAcontent=0.9 - > ${pola_genomicA_filtered_dir}${sample}_dedup_uro_minus_polyAsites.bed --discardInternallyPrimed --genomeFasta='/home/sonas/beegfs/cellranger/ref/refdata-gex-GRCh38-2020-A/fasta/genome.fa' --minUpMisPrimeAlength=10

done


### Now running analysis on 'plus' strand ###

cd ${strand_dir}

f=*plus.bam
for bam_file in ${f[*]};do

sample=${bam_file%_dedup_uro_plus.bam}

samtools view ${bam_file} |/home/sonas/beegfs/apa_tinglab/01_polyAseq/files_modified_for_scAPA/samToPolyA.pl --minClipped=10 --minAcontent=0.9 - > ${pola_genomicA_filtered_dir}${sample}_dedup_uro_plus_polyAsites.bed --discardInternallyPrimed --genomeFasta='/home/sonas/beegfs/cellranger/ref/refdata-gex-GRCh38-2020-A/fasta/genome.fa' --minUpMisPrimeAlength=10

done

#samtools view ${bam_file} |/home/sonas/beegfs/apa_tinglab/01_polyAseq/files_modified_for_scAPA/samToPolyA.pl --minClipped=10 --minAcontent=0.9 - > ${polya_dir}${sample}_dedup_uro_plus_polyAsites.bed




