#!/usr/bin/bash

#SBATCH --mail-type=END
#SBATCH --mail-user=sonas@ccf.org
#SBATCH --job-name=final_peak_ref
#SBATCH -n 20
#SBATCH --mem-per-cpu=3900

#SBATCH -o peak_ref%J.out
#SBATCH -e peak_ref_%J.err


### Load modules ###

module load R/4.1.0

### Initialize variables ###

extn=25
#extn=50
#extn=75
#extn=100
#extn=125
#extn=149


outdir='/home/sonas/beegfs/APA/scPASU/output/2_PeakRef/2d_TU_annotated/'
fprefix='u10_uro_l300iA10mm3_e'${extn}

script_dir='/home/sonas/beegfs/APA/scPASU/scripts/'
${script_dir}update_peak_ref.R ${peak_ref_dir} ${outdir} ${fprefix}



### Archive ###

script_dir='/home/sonas/beegfs/my_scripts/pipeline/my_scAPA/peak_ref/'
bash ${script_dir}create_final_peak_ref.sh ${peak_ref_dir} ${outdir} ${fprefix}


#peak_ref_dir='/home/sonas/beegfs/APA/scAPA/ureter10/assign_tu/u10_uro_unfiltered_3polya/'
#fprefix='u10_uro_unfiltered_3polya'


#peak_ref_dir='/home/sonas/beegfs/APA/scAPA/ureter10/assign_tu/output_3polyA/'
#_A10mm0_3polya'

#peak_ref_dir='/home/sonas/beegfs/APA/scAPA/ureter10/assign_tu/l300i_A6mm0_3polyA/'
#fprefix='u10_uro_final_peak_ref_l300i_A6mm0_3polya'

#peak_ref_dir='/home/sonas/beegfs/APA/scAPA/ureter10/assign_tu/l300i_A10mm3_3polyA/'
#fprefix='u10_uro_final_peak_ref_l300i_A10mm3_3polya'

# Run these
peak_ref_dir='/home/sonas/beegfs/APA/scAPA/ureter10/assign_tu/l250i_A10mm3_3polyA/'
fprefix='u10_uro_final_peak_ref_l250i_A10mm3_3polya'

peak_ref_dir='/home/sonas/beegfs/APA/scAPA/ureter10/assign_tu/l200i_A10mm3_3polyA/'
fprefix='u10_uro_final_peak_ref_l200i_A10mm3_3polya'

#peak_ref_dir='/home/sonas/beegfs/APA/scAPA/ureter10/assign_tu/l300i_A10mm2_3polyA/'
#fprefix='u10_uro_final_peak_ref_l300i_A10mm2_3polya'

# Peak ref after TU assignment
script_dir='/home/sonas/beegfs/my_scripts/pipeline/my_scAPA/peak_ref/'
bash ${script_dir}create_final_peak_ref.sh ${peak_ref_dir} ${outdir} ${fprefix}

# TSS filtered ref - add code for calling respective script here


# Create SAF format ref too
peak_ref_file=${outdir}${fprefix}_updated_tss_filt.txt

Rscript ${script_dir}create_saf.R ${outdir} ${peak_ref_file} ${fprefix}
