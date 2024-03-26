#!/bin/bash

source ~/miniconda3/etc/profile.d/conda.sh
cd /gpfs/fs7/grdi/genarcc/wp3/judsonb/chinook_hdGBS
conda activate ../hdgbs_env

VCF1="populations.snps.vcf.gz"

vcftools --gzvcf "$VCF1" --hardy --depth --missing-indv --missing-site --site-mean-depth --plink

plink --file --indep-pairwise 50 10 0.2

