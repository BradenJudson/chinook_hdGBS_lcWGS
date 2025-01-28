#!/bin/bash


#SBATCH --job-name=otsh_array_impute
#SBATCH --partition=standard
#SBATCH --open-mode=append
#SBATCH --output=03a_imputed_vcf/00_log_folder/beagle_chr.%A_%a.out
#SBATCH --error=03a_imputed_vcf/00_log_folder/beagle_chr.%A_%a.err
#SBATCH --mail-user=Timothy.Healy@dfo-mpo.gc.ca
#SBATCH --time=36:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=110GB
#SBATCH --account=grdi_genarcc
#SBATCH --array=1-34


# Source conda environment.

source ~/bioinfo/mambaforge/etc/profile.d/conda.sh
cd /gpfs/fs7/grdi/genarcc/wp3/healyt/chinook/imputation
conda activate ../../envs/beagle4

# Define global options

REGION=`cat ./chinook_chrs_bp.txt | head -n $SLURM_ARRAY_TASK_ID | tail -n 1`
INPUT="02_input_vcf"
OUTPUT="03a_imputed_vcf"
JAVA_OPTS="-Xmx80G"
NCPU=64

# Impute missing genotypes

name=$(basename $(ls -1 "$INPUT"/*.vcf.gz | perl -pe 's/.vcf.gz//g'))

beagle "$JAVA_OPTS" \
        gl="$INPUT"/"$name".vcf.gz \
        out="$OUTPUT"/"$name"_"$REGION".imputed \
        chrom="$REGION" \
        gprobs=true \
        nthreads=$NCPU
