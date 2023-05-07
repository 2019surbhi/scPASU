#!/usr/bin/bash

#SBATCH --mail-type=END
#SBATCH --mail-user=sonas@ccf.org
#SBATCH --job-name=polya_m
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=3900

#SBATCH -o polyA_m.out
#SBATCH -e polyA_m.err

# conda activate scPASU

export LD_LIBRARY_PATH=/cm/local/apps/gcc/8.2.0/lib64/

module load gcc/8.2.0
module load samtools/1.9

export PATH='/home/sonas/beegfs/miniconda3/envs/scPASU/lib/perl5/site_perl/5.22.0/Bio/DB':$PATH
export PERL5LIB='/home/sonas/beegfs/miniconda3/envs/scPASU/lib/perl5/site_perl/5.22.0/':$PATH

fasta='/home/sonas/beegfs/ref/refdata-gex-GRCh38-2020-A/fasta/genome.fa'

strand='minus'
#strand='plus'

bam_dir='/home/sonas/beegfs/APA/scPASU/output/1_process_bam/1e_merged_bam/'
polya_dir='/home/sonas/beegfs/APA/scPASU/output/2_PeakRef/2b_polya/'

bam=dedup_u10_uro_clean_filtered_${strand}.bam
bed= u10_uro_polyAsites_${strand}.bed
script_dir='/home/sonas/beegfs/APA/scPASU/scripts/'


### Run polyA script ###

cd ${strand_dir}

echo Finding polyA sites in ${bam} and writing output to ${bed}

samtools view ${bam_dir}${bam} | ${script_dir}samToPolyA.pl --minClipped=10 --minAcontent=0.9 - > ${polya_dir}${bed}






### Archive ###

#module load perl/5.34.0

export PATH='/mnt//beegfs/sonas/miniconda3/envs/perl/lib/perl5/site_perl/5.22.0/Bio/DB/':$PATH
export PERL5LIB="/mnt//beegfs/sonas/miniconda3/envs/perl/lib/perl5/site_perl/5.22.0/:"$PATH

fasta='/home/sonas/beegfs/ref/refdata-gex-GRCh38-2020-A/fasta/genome.fa'
script_dir='/home/sonas/beegfs/my_scripts/perl/'

#strand_dir='/home/sonas/beegfs/cellranger/mapped_data/ureter10/ureter10_uro/scAPA_1/temp/dedup_uro_bam_original_name/'
#polya_dir='/home/sonas/beegfs/results/scAPA/my_scAPA/ureter_uro/polya_0.8_merged_BAM/'
#polya_genomicA_filtered_dir='/home/sonas/beegfs/results/scAPA/my_scAPA/ureter_uro/polya_genomicA_filtered_0.8_merged_BAM/'
#script_dir='/home/sonas/beegfs/my_scripts/perl/'

#m='deduped_u10_uro_merged_minus.bam'
#p='deduped_u10_uro_merged_plus.bam'

#####--#####



### Specify INPUTS and OUTPUTS ###

#strand='minus'
strand='plus'

# 1. Finding polyA sites - dedup only BAM

#strand_dir='/home/sonas/beegfs/APA/scAPA/ureter10/macs2_run/macs2_input/dedup_bam/'
#polya_dir='/home/sonas/beegfs/APA/scAPA/ureter10/polya_dir/dedup_only/'
#polya_genomicA_discarded_dir='/home/sonas/beegfs/APA/scAPA/ureter10/polya_dir/dedup_genomicA_discarded/'

#bam=dedup_u10_uro_merged_${strand}.bam

#bed=u10_uro_polyAsites_${strand}.bed
#bed_f=u10_uro_genomicA_discarded_polyAsites_${strand}.bed

# 2. Finding polyA sites - dedup and genomicA filtered BAM

strand_dir='/home/sonas/beegfs/APA/scAPA/ureter10/macs2_run/macs2_input/dedup_genomicA_filtered_bam/'
#polya_dir='/home/sonas/beegfs/APA/scAPA/ureter10/polya_dir/dedup_genomicA_filtered/'
pola_genomicA_discarded_dir='/home/sonas/beegfs/APA/scAPA/ureter10/polya_dir/dedup_genomicA_filtered_genomicA_discarded/'

bam=dedup_u10_uro_genomicA_filtered_merged_${strand}.bam

#bed=u10_uro_genomicA_filtered_polyAsites_${strand}.bed
bed_f=u10_uro_genomicA_filtered_genomicA_discarded_polyAsites_${strand}.bed


### Run polyA script ###

cd ${strand_dir}

### A) polyA sites only ###
#echo Finding polyA sites in ${bam} and writing output to ${bed}

#samtools view ${strand_dir}${bam} | ${script_dir}samToPolyA.pl --minClipped=10 --minAcontent=0.9 - > ${polya_dir}${bed}


### B) polyA sites with genomicA discarded ###
echo Finding polyA sites while discarding genomicA in ${bam} and writing output to ${bed_f}

samtools view ${strand_dir}${bam} | ${script_dir}samToPolyA.pl --minClipped=10 --minAcontent=0.9 - > ${polya_genomicA_discarded_dir}${bed_f} --discardInternallyPrimed --genomeFasta=${fasta} --minUpMisPrimeAlength=10





