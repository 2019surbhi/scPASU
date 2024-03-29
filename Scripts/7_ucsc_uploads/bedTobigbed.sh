#!/usr/bin/bash

#SBATCH --mail-type=END
#SBATCH --mail-user=sonas@ccf.org
#SBATCH --job-name=3polya_bedTobigbed
#SBATCH -n 40
#SBATCH --mem-per-cpu=3900

#SBATCH -o bedTobigbed.out
#SBATCH -e bedTobigbed.err

module load R/4.1.0 
module load python/3.9.5
#module load bedtools/2.29.0


bb_path='/home/sonas/beegfs/tools/ucsc-tools/'
ref_dir='/home/sonas/beegfs/ref/refdata-gex-GRCh38-2020-A/fasta/'


#dir='/home/sonas/beegfs/APA/scAPA/ureter10/peak_ref_final/'
#bed='u10_uro_unfiltered_3polya_updated_new.bed'
#bed_prefix='u10_uro_p1_unfiltered_3polya'

#dir='/home/sonas/beegfs/APA/scAPA/ureter10/rerun/peak_universe_w_extn/final_ref/'
#dir='/home/sonas/beegfs/APA/scPASU/output/3_RefinePeakRef/3b_TSS_filtered/'

#dir='/home/sonas/beegfs/APA/scPASU/output/2_PeakRef/2e_TSS_filtered/'
#bed='u10_uro_l300iA10mm3_e100_tss_filt.bed'
#bed_prefix='u10_uro_l300iA10mm3_e100_tss_filt'

#dir='/home/sonas/beegfs/APA/scPASU/input/'
#bed='pdb_hg38.bed'
#bed_prefix='pdb_hg38'

dir='/home/sonas/beegfs/APA/scPASU/output/2_PeakRef/2e_TSS_filtered/'
bed='u10_uro_l300iA10mm3_e100_updated_tss_filt_final.bed'
bed_prefix='u10_uro_l300iA10mm3_e100_updated_tss_filt_final'


outdir='/home/sonas/beegfs/APA/scPASU/output/7_UCSC_uploads/bigbed/'

# Use R to read the text file, subset relevant columns and then save it using write.table() with sep='\t', col.names=FALSE,row.names=FALSE and quote=FALSE #
# col saved and ordered: chr,start,end,final_annotation, score, strand


cd ${dir}

cat ${bed_prefix}.bed | grep - > ${dir}${bed_prefix}_minus.bed
cat ${bed_prefix}.bed | grep + > ${dir}${bed_prefix}_plus.bed

# Sort

sort -k 1,1 -k2,2n ${dir}${bed_prefix}_minus.bed > ${dir}${bed_prefix}_minus_sorted.bed
sort -k 1,1 -k2,2n ${dir}${bed_prefix}_plus.bed > ${dir}${bed_prefix}_plus_sorted.bed


# Convert to bigBed
bedToBigBed ${dir}${bed_prefix}_minus_sorted.bed ${ref_dir}chr.sizes ${outdir}${bed_prefix}_minus.bb
bedToBigBed ${dir}${bed_prefix}_plus_sorted.bed  ${ref_dir}chr.sizes ${outdir}${bed_prefix}_plus.bb



