#!/bin/bash

#SBATCH --job-name=vcf_conv
#SBATCH --account=grdi_genarcc
#SBATCH --time=24:00:00


source ~/miniconda3/etc/profile.d/conda.sh
cd /gpfs/fs7/grdi/genarcc/wp3/judsonb/chinook_lcWGS/07f_angsd_maf_maf001_nooutliers
conda activate ../../lcwgs_env

CONCAT_VCF="chinook_concat_maf1_m15_nout.vcf.gz"
MAF=0.05

# Ouputs compressed vcf file from 10mb subset bcfs.
bcftools concat *.bcf -O z -o ../09_vcf/"$CONCAT_VCF" --threads 32

cd ../09_vcf

# Subset concatenated vcf to a MAF of 0.05.
bcftools view -O z -i 'MAF[0]>'$MAF'' "$CONCAT_VCF" > chinook_filtered_maf5.vcf.gz


