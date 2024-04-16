#!/bin/bash


#SBATCH --job-name=hdgbs_populations
#SBATCH --partition=standard
#SBATCH --mail-user=Braden.Judson@dfo-mpo.gc.ca
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=20
#SBATCH --comment="registry.maze.science.gc.ca/ssc-hpcs/generic-job:ubuntu22.04"
#SBATCH --account=grdi_genarcc


STACKS_FOLDER="06_stacks_nopcrdedup"
INFO_FILES_FOLDER="01_info_files"
POP_MAP="population_map.txt"
NUM_CPU=20

source ~/.profile

populations -P "$STACKS_FOLDER" \
        -M "$INFO_FILES_FOLDER"/"$POP_MAP" \
        -t "$NUM_CPU" -p 2 -r 0.6 -R 0.6 \
        --ordered-export \
        --write-random-snp \
        --fasta-loci \
        --vcf \
        --plink \
        -O 08_populations_nopcrdedup_rSNP_m60