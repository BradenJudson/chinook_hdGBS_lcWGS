#!/bin/bash


#SBATCH --job-name=hdgbs_fastqc
#SBATCH --partition=standard
#SBATCH --mem-per-cpu=5G
#SBATCH --ntasks=1
#SBATCH --mail-user=Braden.Judson@dfo-mpo.gc.ca
#SBATCH --time=05:00:00
#SBATCH --cpus-per-task=1
#SBATCH --comment="registry.maze.science.gc.ca/ssc-hpcs/generic-job:ubuntu22.04"
#SBATCH --account=grdi_genarcc


# Source conda environment.
source ~/miniconda3/etc/profile.d/conda.sh
cd /gpfs/fs7/grdi/genarcc/wp3/judsonb/chinook_hdGBS
conda activate ../hdgbs_env

# Multi-thread across 6 cores.
fastqc -o 02_reads/trimmed_fastqc/ -t 6 02_reads/trimmed/*fastq.gz
