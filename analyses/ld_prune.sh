#!/bin/bash

#SBATCH --job-name=lcwgs_ld
#SBATCH --partition=standard
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5G
#SBATCH --account=grdi_genarcc

source ~/miniconda3/etc/profile.d/conda.sh
cd /gpfs/fs7/grdi/genarcc/wp3/judsonb/chinook_lcWGS/09_vcf
conda activate ../../lcwgs_env

VCF="chinook_lcwgs_maf1_m15_n361_imputed.vcf.gz"
WINDOWSIZE=50
WINDOWSTEP=10
RSQ=0.4

plink --vcf "$VCF" \
        --double-id \
        --aec \
        --set-missing-var-ids @:# \
        --indep-pairwise $WINDOWSIZE $WINDOWSTEP $RSQ \
        --out imputed_vcf_ld

perl -pe 's/:/\t/g' imputed_vcf_ld.prune.in > imputed_vcf_prune_keep

vcftools --gzvcf "$VCF" \
        --positions imputed_vcf_prune_keep \
        --recode --stdout | bgzip - > chinook_lcwgs_maf1_m15_n361_imputed_pruned.vcf.gz
