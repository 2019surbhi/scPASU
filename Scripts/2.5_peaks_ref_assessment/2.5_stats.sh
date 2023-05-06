#!/usr/bin/bash

#SBATCH --mail-type=END
#SBATCH --mail-user=sonas@ccf.org
#SBATCH --job-name=stats
#SBATCH -c 20
#SBATCH --mem 70G

#SBATCH -o stats_m.out
#SBATCH -e stats_m.err

module load R

strand='plus'
#strand='minus'

min_polya=3

extn=25
#extn=50
#extn=75
#extn=100
#extn=125
#extn=149

fprefix=u10_uro_l300iA10mm3_e${extn}

script_dir='/home/sonas/beegfs/APA/scPASU/scripts/'
dir='/home/sonas/beegfs/APA/scPASU/output/2_PeakRef/2c_intersect/e'${extn}/${strand}/

peaks_dir=${dir}peaks_dir/
polya_dir=${dir}polya_dir/
int_dir=${dir}int_dir/
pr_dir=${dir}pr_dir/
pr_table_dir=${dir}pr_table_dir/
peaks_count_dir=${dir}peaks_count_dir/
polya_count_dir=${dir}polya_count_dir/
final_filtered_ref_table=${dir}final_filtered_ref_table/

cd ${dir}
mkdir stats

stats=${dir}stats/
outdir='/home/sonas/beegfs/APA/scPASU/output/2.5_PeakRefAssessment/stats/'

# Get summary count for all intermediate steps

# total peaks
cd ${peaks_dir}
f=*sorted*
wc -l $f > ${stats}total_peaks_stats_${strand}.txt
 
# peaks supported by polya

cd ${peaks_count_dir}
f=*
wc -l $f > ${stats}peaks_supported_by_polya_${strand}.txt

# peaks supported my min_polya reads
cd ${peaks_count_dir}
f=*
for file in ${f[*]};do
count=$(cat ${file} | awk '{if($1>=3) {print $1"\t"$2};}' | wc -l)
printf '%s\t%s\n' ${count} ${file}>> ${stats}peaks_supported_by_min_polya_${strand}.txt
done

# peaks not supported by polya reads 
cd ${int_dir}
f=*tagged_sorted_${strand}.bed
for files in ${f[*]};do count=$(cat $files | grep non_polya | wc -l); printf "%s\t%s\n" ${count} ${files} >> ${stats}peaks_not_supported_by_polya_${strand}.txt ; done

# peaks not supported by min polya reads 

cd ${peaks_count_dir}
f=*
for file in ${f[*]};do
count=$(cat ${file} | awk '{if($1<3) {print $1"\t"$2};}' | wc -l)
printf '%s\t%s\n' ${count} ${file}>> ${stats}peaks_not_supported_by_min_polya_${strand}.txt
done


#total polya reads
cd ${polya_dir}
f=*sorted*
wc -l $f > ${stats}total_polya_stats_${strand}.txt
 
#  polya reads supporting peaks
cd ${int_dir}
f=*tagged_sorted_${strand}.bed
for files in ${f[*]};do 
count=$(cat $files | grep -w polya | wc -l)
printf "%s\t%s\n" ${count} ${files}  >> ${stats}polya_reads_supporting_peaks_${strand}.txt
done


#  min polya reads supporting peaks
cd ${polya_count_dir}
f=*count_sorted_${strand}.txt
for files in ${f[*]};do 
count=$(cat ${file} | awk '{if($1>3) {print $1"\t"$2};}' | wc -l)
printf "%s\t%s\n" ${count} ${files}  >> ${stats}min_polya_reads_supporting_peaks_${strand}.txt
done


Rscript ${script_dir}2.5_combine_stats.R ${stats} ${outdir} ${strand}




#### Archive #####

#dir='/home/sonas/beegfs/APA/scAPA/ureter10/u10_uro_peak_universe/u10_uro_peak_universe_minus/'
#dir=/home/sonas/beegfs/APA/scAPA/ureter10/peak_universe/peak_universe_dedup_only/${strand}/
#dir=/home/sonas/beegfs/APA/scAPA/ureter10/peak_universe/peak_universe_l300_A10mm3/${strand}/
#dir=/home/sonas/beegfs/APA/scAPA/ureter10/peak_universe/peak_universe_l300_i_A6mm0/${strand}/
#dir=/home/sonas/beegfs/APA/scAPA/ureter10/peak_universe/peak_universe_l300_i_A10mm3/${strand}/
#dir=/home/sonas/beegfs/APA/scAPA/ureter10/peak_universe/peak_universe_l250_i_A10mm3/${strand}/
#dir=/home/sonas/beegfs/APA/scAPA/ureter10/peak_universe/peak_universe_l200_i_A10mm3/${strand}/
#dir=/home/sonas/beegfs/APA/scAPA/ureter10/peak_universe/peak_universe_l300_i_A10mm2/${strand}/

