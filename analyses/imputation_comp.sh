#!/bin/bash

#SBATCH --account=grdi_genarcc
#SBATCH --partition=standard
#SBATCH --job-name=imp_vcf
#SBATCH --time=24:00:00
#SBATCH --array=1-366
#SBATCH --mem=60GB

source ~/miniconda3/etc/profile.d/conda.sh
cd /gpfs/fs7/grdi/genarcc/wp3/judsonb/chinook_lcWGS/12b_imputation_full
conda activate ../../lcwgs_env

IMPVCF="../../../healyt/chinook/imputation/04_merged_vcf/chinook_concat_maf5_m15_nout.imputed.vcf.gz"
ORIVCF="../../../healyt/chinook/imputation/02_input_vcf/chinook_concat_maf5_m15_nout.vcf.gz"
HDGBS="../../chinook_hdGBS/09c_filter_branch2/snps_maf001.vcf.gz"

#cp "$IMPVCF" ./chinook_imputed.vcf.gz
#cp "$ORIVCF" ./chinook_original.vcf.gz
#cp "$HDGBS"  ./hdgbs_snps.vcf.gz

HDGBS_file="hdgbs_snps.vcf.gz"
Imp_VCF="chinook_imputed.vcf.gz"
Ori_VCF="chinook_original.vcf.gz"

#tabix -p vcf "$HDGBS_file"
#tabix -p vcf "$Imp_VCF"
#tabix -p vcf "$Ori_VCF"

# Bam list and paths for each sample. These are the sample names in the lcWGS datasets.
bam=`cat shared_samples.txt | head -n $SLURM_ARRAY_TASK_ID | tail -n 1`

# Individual sample names. These are the names of the samples in the hdGBS dataset, but are also used in naming.
ind=`cat shared_samples.txt | head -n $SLURM_ARRAY_TASK_ID | tail -n 1 | perl -pe 's/^.*i5_[0-9]*\.//g' | perl -pe 's/\.dedup.*$//g' -`

        # For the original (non-imputed) lcWGS dataset, isolate the missing SNPs on an individual basis.
        awk 'BEGIN {OFS="\t"} {if ($6 == 1) {print $1, $2}}' <(vcftools --gzvcf "$Ori_VCF" --indv "$bam" --missing-site --stdout) > missing_snps/"$ind"_missing_snps.txt

        # Filter the hdGBS dataset to retain only those SNPs that were originally missing in the non-imputed lcWGS VCF and retain genotypic information.
        awk '{if ($5 != "\.\/\.") print $0}' <(bcftools query -s "$ind" -R missing_snps/"$ind"_missing_snps.txt -f '%CHROM\t%POS\t%REF\t%ALT\t[ %GT]\n' hdgbs_snps.vcf.gz) > hdgbs_geno/"$ind"_hdgbs_geno.txt

        # Subset the imputed lcWGS VCF to have the same SNPs as present in the individual hdGBS files and output genotypic information for each individual.
        bcftools query -s "$bam" -R <(awk 'BEGIN {OFS="\t"} {print $1, $2}' hdgbs_geno/"$ind"_hdgbs_geno.txt) -f '%CHROM\t%POS\t%REF\t%ALT\t[ %GT]\n' "$Imp_VCF" > imp_geno/"$ind"_imp_geno.txt

# Isolate SNPs shared between the hdGBS and imputed lcWGS datasets.
gawk 'NR==FNR {a[$1,$2]++; next} (($1,$2) in a) {print $0}' hdgbs_geno/"$ind"_hdgbs_geno.txt imp_geno/"$ind"_imp_geno.txt > comp_files/"$ind"_shared.txt

# Count the SNPs per file and save them in a text file.
wc -l comp_files/* > shared_snps_imp_hdgbs.txt

# Extract individual-level allele frequencies for the SNPs identified above.
vcftools --gzvcf "$HDGBS_file" --indv "$ind" --positions comp_files/"$ind"_shared.txt --freq --out hdgbs_freq/"$ind"_freq_GLs.txt
vcftools --gzvcf "$Imp_VCF"    --indv "$bam" --positions comp_files/"$ind"_shared.txt --freq --out imp_freq/"$ind"_freq_GLs.txt

# Extract global allele frequencies for imputed and non-imputed low coverage datasets.
vcftools --gzvcf "$Imp_VCF" --freq --out imp_freq/global_imputed_freqGLs.txt
vcftools --gzvcf "$Ori_VCF" --freq --out imp_freq/global_original_freqGLs.txt