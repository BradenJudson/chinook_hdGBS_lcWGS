#!/bin/bash


source ~/miniconda3/etc/profile.d/conda.sh
cd /gpfs/fs7/grdi/genarcc/wp3/judsonb/chinook_hdGBS/09_filtering
conda activate ../../hdgbs_env

cp ../08_populations_nopcrdedup_rSNP_m60/populations.snps.vcf populations.snps.vcf

python extract_snp_with_max_maf.py "$FDIR"/unfiltered_populations.snps.vcf "$FDIR"/"$VCF1"

VCF1="populations.snps.vcf"

# Remove low quality genotypes.
vcftools --vcf "$VCF1" --minGQ 20 --min-alleles 2 --max-alleles 2 --remove-indels --recode --stdout | gzip -c > highqsnps_m60.vcf.gz
bcftools stats highqsnps_m60.vcf.gz > highqsnps_m60.txt

# Remove loci that aren't biallelic.
vcftools --vcf highqsnps_m60.vcf --remove-indels --min-alleles 2 --max-alleles 2 --recode --stdout > snps_biallelic_m60.vcf
bcftools stats snps_biallelic_m60.vcf > snps_biallelic_m60.txt

# Remove individuals with too few reads.
vcftools --gzvcf highqsnps_m60.vcf.gz --remove ind_fails.txt --recode --stdout | gzip -c > filtered_indvs_m60.vcf.gz
bcftools stats filtered_indvs_m60.vcf.gz > filtered_indvs_m60.txt

# Remove loci with coverage outside of thresholds.
vcftools --gzvcf filtered_indvs_m60.vcf.gz --min-meanDP 3 --max-meanDP 100 --recode --stdout | gzip -c > snps_coverage_m60.vcf.gz
bcftools stats snps_coverage_m60.vcf.gz > snps_coverage_m60.txt

# Remove loci not present in 70% of individuals globally.
vcftools --gzvcf snps_coverage_m60.vcf.gz --max-missing 0.3  --recode --stdout | gzip -c > max70missing_m60.vcf.gz
bcftools stats max70missing_m60.vcf.gz > max70missing_m60.txt

# Remove loci with MAF < 5%.
vcftools --gzvcf max70missing_m60.vcf.gz --maf 0.01 --recode --stdout | gzip -c > snps_maf1_m60.vcf.gz
bcftools stats snps_maf1_m60.vcf.gz > snps_maf1_m60.txt

# Next use HDPLot in R for removing paralogs and implementing a MAF cutoff of 5%.