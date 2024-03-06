#!/bin/bash


#SBATCH --job-name=hdgbs_cutadapt
#SBATCH --partition=standard
#SBATCH --mem-per-cpu=5G
#SBATCH --ntasks=1
#SBATCH --mail-user=Braden.Judson@dfo-mpo.gc.ca
#SBATCH --time=16:00:00
#SBATCH --cpus-per-task=64
#SBATCH --comment="registry.maze.science.gc.ca/ssc-hpcs/generic-job:ubuntu22.04"
#SBATCH --account=grdi_genarcc


# Source conda environment.
source ~/miniconda3/etc/profile.d/conda.sh
cd /gpfs/fs7/grdi/genarcc/wp3/judsonb/chinook_hdGBS
conda activate ../hdgbs_env

NCPU=64

# Load GNU parallel.
. ssmuse-sh -x main/opt/parallel/parallel-202109022

# Export for parallelization.
export cutadapt

ls -1 02_reads/*_R1.fastq.gz |
         sed 's/_R1\.fastq\.gz//g' |
         sed 's/02\_reads\///g' |
parallel -j $NCPU cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
         -g AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
         -o 02_reads/trimmed/{}"_R1.fastq.gz" \
         -p 02_reads/trimmed/{}"_R2.fastq.gz" \
         -e 0.2 \
         -m 50 \
            02_reads/{}_R1.fastq.gz \
            02_reads/{}_R2.fastq.gz





