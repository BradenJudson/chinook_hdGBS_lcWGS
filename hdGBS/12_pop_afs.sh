#!/bin/bash

source ~/miniconda3/etc/profile.d/conda.sh
cd /gpfs/fs7/grdi/genarcc/wp3/judsonb/chinook_hdGBS/14_pop_frq
conda activate ../../hdgbs_env

SAMPLEINFO="../01_info_files/hdGBS_sampleinfo_fixed.csv"
VCF="../13_imputation/snps_maf005_imputed.vcf.gz"


POPLIST=`awk -F "," '{print $3}' ../01_info_files/hdGBS_sampleinfo_fixed.csv | sort | uniq`

for POP in $POPLIST
do
        echo $POP

        awk -v pop="$POP" -F"," '{ if ($3 == pop) {print $4}}' "$SAMPLEINFO" > "$POP"_indvs.txt

        vcftools --gzvcf "$VCF" --keep "$POP"_indvs.txt --freq --out "$POP"_freqs

done
