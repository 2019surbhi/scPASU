#!/usr/bin/bash

## This script creates peak ref using MACS2 and polya reads outputs ##

# Author: Surbhi Sona

## User inputs ##

OPTIND=1

while getopts "m:p:o:s:c:n:e:t:" opt;do
  case $opt in
  m)
    macs_dir=$OPTARG
    ;; #-m path to MACS2 peaks output
  p)
    polya_reads_dir=$OPTARG
    ;; #-p path to polya dir
  o)
    out_dir=$OPTARG
    ;; #-o specifies output directory for peak ref
  s) 
    strand=$OPTARG
    ;; #-s specify strand as plus or minus
  c)
    chr=$OPTARG
    ;; # -c specify chr (1-22 or X or Y)
  n)
    min_polya=$OPTARG
    ;; #-n minimum polya reads that a peak must be supported by
   \?) echo "Option ${opt} not valid";;   
  e)
    extn=$OPTARG\
    ;; #-e extension length (3' of peak extended by this length) 
  t)
    script_dir=$OPTARG
    ;; #-t script dir path 
  esac
done

### Set paths ###


cd ${out_dir}
outdir=${out_dir}${strand}/
mkdir ${outdir}

peaks_dir=${outdir}peaks_dir/
polya_dir=${outdir}polya_dir/
int_dir=${outdir}int_dir/
pr_dir=${outdir}pr_dir/
pr_table_dir=${outdir}pr_table_dir/
peaks_count_dir=${outdir}peaks_count_dir/
polya_count_dir=${outdir}polya_count_dir/
ref_table=${outdir}ref_table/
filtered_ref_table=${outdir}filtered_ref_table/

### 0. Prepare files ###

#(a) create dir

cd ${outdir}

mkdir ${peaks_dir}
mkdir ${polya_dir}
mkdir ${int_dir}
mkdir ${pr_dir}
mkdir ${pr_table_dir}
mkdir ${peaks_count_dir}
mkdir ${polya_count_dir}
mkdir ${ref_table}
mkdir ${filtered_ref_table}
 

# Create chr file (currently doing it manually to exclude non standard chr)

#awk '{print $1}' ${peaks_file} | sort | uniq | sort -V -k1,1 > ${dir}chr_${strand}.txt

# (b) Specify peaks and polya  files

cd ${macs_dir}
peaks_file=${macs_dir}*_${strand}_peaks.narrowPeak
polya_file=${polya_reads_dir}*polyAsites*_${strand}.bed


# (c) Split peak columns by chr, subset desired col and add strand column and create create x base extension at 3' end for peaks , where x is defined by user as $extn
# I intend to try 100, 125 and 149

echo creating peak and polya files by chr

if [ ${strand} == 'plus' ]
then
   cat ${peaks_file} | grep -w ${chr} | awk 'OFS="\t" {print $1,$2,$3=$3+ '$extn',$4,$5="+"}' > ${peaks_dir}peaks_${chr}_${strand}.bed
else
   cat ${peaks_file} | grep -w ${chr} | awk 'OFS="\t" {print $1,$2=$2- '$extn',$3,$4,$5="-"}' > ${peaks_dir}peaks_${chr}_${strand}.bed
fi
 
# subset poly table too

cat ${polya_file} | grep -w ${chr} | awk 'OFS="\t" {print $1,$2,$3,$4,$5,$6}' > ${polya_dir}polya_${chr}_${strand}.bed


# (e) sort 

echo sort 

cat  ${peaks_dir}peaks_${chr}_${strand}.bed | sort -k2,2n -k3,3n >  ${peaks_dir}peaks_${chr}_sorted_${strand}.bed
cat  ${polya_dir}polya_${chr}_${strand}.bed | sort -k2,2n -k3,3n >  ${polya_dir}polya_${chr}_sorted_${strand}.bed

### 1a. Run intersect ###

echo intersect

bedtools intersect -a ${peaks_dir}peaks_${chr}_sorted_${strand}.bed -b ${polya_dir}polya_${chr}_sorted_${strand}.bed -wa -wb -loj > ${int_dir}peaks_int_polya_${chr}_${strand}.bed


### 1b. Sort intersect table ###
cat ${int_dir}peaks_int_polya_${chr}_${strand}.bed | sort -k6,6n  > ${int_dir}peaks_int_polya_${chr}_sorted_${strand}.bed

### 2. Process table ###

### 2a. Tag non-poly peaks

awk  'FS=OFS="\t" {if($6==".") {$12="non_polya"} else {$12="polya"}; {print $0}}' ${int_dir}peaks_int_polya_${chr}_sorted_${strand}.bed > ${int_dir}peaks_int_polya_${chr}_tagged_${strand}.bed

### 2b. Sort by polya column

cat ${int_dir}peaks_int_polya_${chr}_tagged_${strand}.bed | sort -k12,12 > ${int_dir}peaks_int_polya_${chr}_tagged_sorted_${strand}.bed

### 3a. Create peaks count
cat ${int_dir}peaks_int_polya_${chr}_tagged_${strand}.bed | grep -w polya | awk 'FS=OFS="\t" {print $4}' | sort | uniq -c | sort -k1,1n > ${peaks_count_dir}peaks_${chr}_count_sorted_${strand}.txt

### 3.b. Create polya read count per peak
cd ${int_dir}
file=${int_dir}peaks_int_polya_${chr}_tagged_${strand}.bed
peaks_file=${peaks_count_dir}peaks_${chr}_count_sorted_${strand}.txt
peaks=$(cat ${peaks_file} | awk '{print $2}')

for peak in ${peaks[*]};do
count=$(cat ${file} | grep -w ${peak} | wc -l)
printf '%s\t%s\n' ${count} ${peak} >> ${polya_count_dir}polya_${chr}_count_${strand}.txt
done

# sort
cd ${polya_count_dir}
cat polya_${chr}_count_${strand}.txt | sort -k1,1n > polya_${chr}_count_sorted_${strand}.txt 

### 3c. Create PR table ###

cat ${int_dir}peaks_int_polya_${chr}_tagged_sorted_${strand}.bed | sort -u -k4,4  > ${pr_table_dir}pr_table_${chr}_${strand}.bed


### 4. Get inferred PR ###

f=$(awk '{print $2}' ${peaks_count_dir}peaks_${chr}_count_sorted_${strand}.txt)

for peaks in ${f[*]};do
cat ${int_dir}peaks_int_polya_${chr}_sorted_${strand}.bed | grep -w ${peaks}  | awk '{if(min==""){min=max=$7}; if($8>max) {max=$8}; if($7< min) {min=$7};} END {print $4"\t"min"\t"max}' >> ${pr_dir}pr_${chr}_${strand}.bed
done

# Sort

cat ${pr_dir}pr_${chr}_${strand}.bed | sort -k2,2n -k3,3n > ${pr_dir}pr_${chr}_sorted_${strand}.bed

### 5. Final PR table ###

Rscript ${script_dir}1a_create_peak_ref.R ${outdir} ${strand} ${chr}

# Filter peaks by min polya reads

bed_file=${ref_table}ref_table_updated_${chr}_${strand}.bed

Rscript ${script_dir}1b_filter_peak_ref.R ${bed_file} ${filtered_ref_table} ${min_polya}

