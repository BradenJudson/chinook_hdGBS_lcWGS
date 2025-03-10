#!/bin/bash


#SBATCH --job-name=lcwgs_fastp
#SBATCH --partition=standard
#SBATCH --mem-per-cpu=5G
#SBATCH --ntasks=1
#SBATCH --mail-user=Braden.Judson@dfo-mpo.gc.ca
#SBATCH --time=16:00:00
#SBATCH --cpus-per-task=16
#SBATCH --comment="registry.maze.science.gc.ca/ssc-hpcs/generic-job:ubuntu22.04"
#SBATCH --account=grdi_genarcc


# Source conda environment.
source ~/miniconda3/etc/profile.d/conda.sh
cd /gpfs/fs7/grdi/genarcc/wp3/judsonb/chinook_lcWGS
conda activate ../lcwgs_env

LENGTH=100
QUAL=20
INPUT="02_reads"
OUTPUT="03_trimmed_reads"
METRICS="11_metrics"
NCPU=16
#SAMPLEFILE=$1

# Load GNU parallel.
. ssmuse-sh -x main/opt/parallel/parallel-202109022

# Export for parallelization.
export fastp

ls "$INPUT"/*_R1.fastq.gz | perl -pe 's/R[12]\.fastq\.gz//g' | \
        parallel -j "$NCPU" \
                fastp -i {}R1.fastq.gz -I {}R2.fastq.gz \
                -o "$OUTPUT"/{/}R1.fastq.gz \
                -O "$OUTPUT"/{/}R2.fastq.gz \
                --length_required="$LENGTH" \
                --qualified_quality_phred="$QUAL" \
                --correction \
                --trim_tail1=1 \
                --trim_tail2=1 \
                --json 11_metrics/{/}.json \
                --html 11_metrics/{/}.html \
                --report_title={/}.html