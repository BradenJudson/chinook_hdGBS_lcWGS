#!/bin/bash

#SBATCH --job-name=gls
#SBATCH --account=grdi_genarcc
#SBATCH --time=48:00:00
#SBATCH --mem=60GB

source ~/miniconda3/etc/profile.d/conda.sh
cd /gpfs/fs7/grdi/genarcc/wp3/judsonb/chinook_lcWGS/12b_imputation_full
conda activate ../../lcwgs_env

IMPVCF="../../../healyt/chinook/imputation/04_merged_vcf/chinook_concat_maf5_m15_nout.imputed.vcf.gz"
ORIVCF="../../../healyt/chinook/imputation/02_input_vcf/chinook_concat_maf5_m15_nout.vcf.gz"

Imp_VCF="chinook_imputed.vcf.gz"
Ori_VCF="chinook_original.vcf.gz"

# Extract genotype probabilities for each genomic position for both datasets.
bcftools query -f "%CHROM\t%POS[\t%GP]\n" "$Imp_VCF" > imputed_GLs.txt
bcftools query -f "%CHROM\t%POS[\t%GP]\n" "$Ori_VCF" > original_GLs.txt

# Determine number of individuals to consider.
cols=`head -1 original_GLs.txt | tr " " '\n' | wc -l`

# Isolate chromosomes and genomic positions per SNP.
cat <(awk '{print $1,$2}' imputed_GLs.txt) > genomic_positions.txt

# Print the highest GP for each locus for each individual across both datasets.
# Formatting in seq is necessary to keep output order consistent with order of individuals.
# Imputed data first:
for i in $(seq -f '%03g' 3 "$cols");
do
        cat <(awk -F "\t" '{ print $'$i' }' imputed_GLs.txt | tr "," "\t" | perl -lane '@a=sort @F;print join "\t", @a' | awk '{ print $3 }') > ./temp_imp/tempGLs_$i.txt
done

# Original, non-imputed data second:
for i in $(seq -f '%03g' 3 "$cols");
do
        cat <(awk -F "\t" '{ print $'$i' }' imputed_GLs.txt | tr "," "\t" | perl -lane '@a=sort @F;print join "\t", @a' | awk '{ print $3 }') > ./temp_ori/tempGLs_$i.txt
done

# Re-add the genomic position data to the outputs of the for-loops above and convert to long form.
paste ./temp_imp/tempGLs_*.txt | awk '{show=0; for (i=1; i<=NF; i++) {if ($i!=0) show=1; col[i]+=$i;}}
                show==1{tr++; for (i=1; i<=NF; i++) vals[tr,i]=$i; tc=NF} END{for(i=1; i<=tr; i++) { for (j=1; j<=tc; j++)
                         { if (col[j]>0) printf("%s%s", vals[i,j], OFS)} print ""; } }' > maximum_imputed_GLs.txt

# Convert to a long-form csv with a header row.
awk '{for(i=3;i<=NF;i++) print $1,$2,$i}' <(paste genomic_positions.txt maximum_imputed_GLs.txt) | sed '1i CHR POS GP' | sed 's/ \+/,/g' > maximum_imputed_GLs_LF.csv

# Re-add the genomic position data to the outputs of the for-loops above and convert to long form.
paste ./temp_ori/tempGLs_*.txt | awk '{show=0; for (i=1; i<=NF; i++) {if ($i!=0) show=1; col[i]+=$i;}}
                show==1{tr++; for (i=1; i<=NF; i++) vals[tr,i]=$i; tc=NF} END{for(i=1; i<=tr; i++) { for (j=1; j<=tc; j++)
                         { if (col[j]>0) printf("%s%s", vals[i,j], OFS)} print ""; } }' > maximum_original_GLs.txt

awk '{for(i=3;i<=NF;i++) print $1,$2,$i}' <(paste genomic_positions.txt maximum_original_GLs.txt) | sed '1i CHR POS GP' | sed 's/ \+/,/g' > maximum_original_GLs_LF.csv