#!/bin/bash

#SBATCH --job-name=lcwgs_ld
#SBATCH --partition=standard
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5G
#SBATCH --account=grdi_genarcc

source ~/miniconda3/etc/profile.d/conda.sh
cd /gpfs/fs7/grdi/genarcc/wp3/judsonb/chinook_lcWGS/09_vcf
conda activate ../../../healyt/envs/vcftools

VCF_original="chinook_lcwgs_maf1_m15_n361_imputed.vcf.gz"
VCF_ldprune="lcwgs_ldpruned_maf005"
BEAGLE_original="../07h_angsd_n361_v2/angsd_c_filtered7M.beagle.gz"
BEAGLE_ldprune="chinook_lcwgs_maf5_ldpruned"
BEAGLE_sites="ld_pruned_snps_beagle_format"
WINDOWSIZE=50
WINDOWSTEP=10
RSQ=0.4

# Identify SNPs in LD and apply consistent position-based formatting.
plink --vcf "$VCF_original" \
        --double-id \
        --aec \
        --set-missing-var-ids @:# \
        --indep-pairwise $WINDOWSIZE $WINDOWSTEP $RSQ \
        --out imputed_vcf_ld

# Modify file containing LD SNPs so vcftools can use it for filtering.
perl -pe 's/:/\t/g' imputed_vcf_ld.prune.in > imputed_vcf_prune_keep.txt

# Prune original VCF to only the sites selected above.
vcftools --gzvcf "$VCF_original" \
        --positions imputed_vcf_prune_keep.txt \
        --recode --stdout | bgzip - > "$VCF_ldprune".vcf.gz

# To trim the complementary genotype likelihood file requires slightly different formatting
# as the first column is CHROMOSOME_POSITION, so we make that file here.
sed 's/\:/_/' imputed_vcf_ld.prune.in > "$BEAGLE_sites".txt

# Read in the original beagle file and subset it down to the sites returned from LD-pruning step above.
cat <(zcat "$BEAGLE_original" | head -n1) <(zcat "$BEAGLE_original" | awk 'NR==FNR{a[$1]=$1OFS$2;next}{$1=a[$1];print}' OFS='\t' "$BEAGLE_sites".txt - | grep -Pv "^\t") | gzip >  "$BEAGLE_ldprune".gz

