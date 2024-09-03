#!/bin/bash
#SBATCH --array=1-453

source ~/miniconda3/etc/profile.d/conda.sh
cd /gpfs/fs7/grdi/genarcc/wp3/judsonb/chinook_lcWGS/12b_imputation_full

IMPVCF="../../../healyt/chinook/imputation/04_merged_vcf/chinook_concat_maf5_m15_nout.imputed.vcf.gz"
ORIVCF="../../../healyt/chinook/imputation/02_input_vcf/chinook_concat_maf5_m15_nout.vcf.gz"

cp "$IMPVCF" "chinook_imputed.vcf.gz"
cp "ORIVCF" "chinook_original.vcf.gz"

Imp_VCF="chinook_imputed.vcf.gz"
Ori_VCF="chinook_original.vcf.gz"

tabix "$Imp_VCF"
tabix "$Ori_VCF"

# Extract genotype probabilities for each genomic position for both datasets.
bcftools query -f "%CHROM\t%POS[\t%GP]\n" "$Imp_VCF" > imputed_GLs.txt
bcftools query -f "%CHROM\t%POS[\t%GP]\n" "$Ori_VCF" > original_GLs.txt

# Isolate chromosomes and genomic positions per SNP to add later on.
cat <(awk '{print $1,$2}' imputed_GLs.txt) > genomic_positions.txt

# Extract maximum genotype probability per SNP. Write individual by individual.
cat <(awk -F "\t" '{ print $'$SLURM_ARRAY_TASK_ID' }' imputed_GLs.txt | tr "," "\t" | perl -lane '@a=sort @F;print join "\t", @a' | awk '{ print $3 }') > ./temp_array_imp/tempGLs_$SLURM_ARRAY_TASK_ID.txt

# Same as above, but requires extra formatting to properly handle numbers in scientific notation.
cat <(awk -F "\t" '{ print $'$SLURM_ARRAY_TASK_ID' }' original_GLs.txt | tr "," "\t" | awk '{for (i=1;i<=NF;i++) if ($i+0 == $i && $i ~ /e/) $i = sprintf("%.10f", $i) } 1' |\
                awk ' {split( $0, a, " " ); asort( a ); for( i = 1; i <= length(a); i++ ) printf( "%s ", a[i] ); printf( "\n" ); }' | awk '{ print $3 }') > ./temp_array_ori/tempGLs_$SLURM_ARRAY_TASK_ID.txt

# Define header: Chromosome name, position, and individual ID. Convert from txt to csv.
HEADER=`cat ../01_info_files/gp_header.txt | sed 's/ \+/,/g'`

# Combine outputs and write to csv.
paste genomic_positions.txt <(paste ./temp_array_imp/*.txt) | sed 's/ \+/,/g' > maximum_imputed_GLs.csv
paste genomic_positions.txt <(paste ./temp_array_ori/*.txt) | sed 's/ \+/,/g' > maximum_original_GLs.csv

# Add header to existing csvs.
sed -i "1i $HEADER" maximum_imputed_GLs.csv
sed -i "1i $HEADER" maximum_original_GLs.csv
