#!/bin/bash

#SBATCH --job-name=geno_probs
#SBATCH --account=grdi_genarcc
#SBATCH --time=4:00:00
#SBATCH --mem=60GB
#SBATCH --array=1-385

source ~/miniconda3/etc/profile.d/conda.sh
cd /gpfs/fs7/grdi/genarcc/wp3/judsonb/chinook_lcWGS/12_imputation_comparison_n385
conda activate ../../lcwgs_env

Imp_VCF="../09_vcf/chinook_lcwgs_maf005_n385_imputed_filtered.vcf.gz"
Ori_VCF="../09_vcf/chinook_lcwgs_maf005_n385.vcf.gz"

# Extract genotype probabilities for each genomic position for both datasets and remove leading tab.
bcftools query -f "[\t%GP]\n" "$Imp_VCF" | sed "s/^[\t]*//" > imputed_GLs.txt
bcftools query -f "[\t%GP]\n" "$Ori_VCF" | sed "s/^[\t]*//" > original_GLs.txt
bcftools query -f "%CHROM\t%POS\n" "$Imp_VCF" > genomic_positions.txt

# Default IDs are the individual bam file names, so we isolate the sample names here (different order for the imputed/non-imputed datasets).
id1=`bcftools query -l ../09_vcf/$Imp_VCF | head -n $SLURM_ARRAY_TASK_ID | tail -n 1 | perl -pe 's/.dedup.clip.bam//g' | perl -pe 's/.dedup.clip.downsampled.bam//g' | sed 's/.*\.//'`
id2=`bcftools query -l ../09_vcf/$Ori_VCF | head -n $SLURM_ARRAY_TASK_ID | tail -n 1 | perl -pe 's/.dedup.clip.bam//g' | perl -pe 's/.dedup.clip.downsampled.bam//g' | sed 's/.*\.//'`

# Isolate chromosomes and genomic positions per SNP to add later on.
cat <(awk '{print $1,$2}' imputed_GLs.txt) > genomic_positions.txt

cat <(awk -F "\t" '{ print $'$SLURM_ARRAY_TASK_ID' }' imputed_GLs.txt | tr "," "\t" | perl -lane '@a=sort @F;print join "\t", @a' | awk '{ print $3 }') > ./temp_array_imp/$id1.txt

cat <(awk -F "\t" '{ print $'$SLURM_ARRAY_TASK_ID' }' original_GLs.txt | tr "," "\t" | awk '{for (i=1;i<=NF;i++) if ($i+0 == $i && $i ~ /e/) $i = sprintf("%.10f", $i) } 1' |\
                awk ' {split( $0, a, " " ); asort( a ); for( i = 1; i <= length(a); i++ ) printf( "%s ", a[i] ); printf( "\n" ); }' | awk '{ print $3 }') | awk '{printf "%.2f\n",$1}' > ./temp_array_ori/$id2.txt

## Combine outputs and write to csv.
paste genomic_positions.txt <(paste ./temp_array_imp/*.txt) | awk '{print NR, $0}' | sed 's/[[:space:]]\+/,/g' > maximum_imputed_GLs.csv
paste genomic_positions.txt <(paste ./temp_array_ori/*.txt) | awk '{print NR, $0}' | sed 's/[[:space:]]\+/,/g' > maximum_original_GLs.csv

# Define header: Row number, chromosome name, position, and individual ID. One for each VCF.
HEADER_IMP=`echo -e "row\nchr\npos\n$(bcftools query -l "$Imp_VCF")" | tr '\n' ',' | sed '$ s/.$//'`
HEADER_ORI=`echo -e "row\nchr\npos\n$(bcftools query -l "$Ori_VCF")" | tr '\n' ',' | sed '$ s/.$//'`

# Add header to existing csvs.
sed -i "1i $HEADER_IMP" maximum_imputed_GLs.csv
sed -i "1i $HEADER_ORI" maximum_original_GLs.csv