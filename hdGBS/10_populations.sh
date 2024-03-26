#!/bin/bash


#SBATCH --job-name=hdgbs_populations
#SBATCH --partition=standard
#SBATCH --mem=100GB
#SBATCH --ntasks=1
#SBATCH --mail-user=Braden.Judson@dfo-mpo.gc.ca
#SBATCH --time=36:00:00
#SBATCH --cpus-per-task=10
#SBATCH --comment="registry.maze.science.gc.ca/ssc-hpcs/generic-job:ubuntu22.04"
#SBATCH --account=grdi_genarcc


STACKS_FOLDER="06_stacks"
INFO_FILES_FOLDER="01_info_files"
POP_MAP="population_map.txt"
NUM_CPU=10

source ~/miniconda3/etc/profile.d/conda.sh
cd /gpfs/fs7/grdi/genarcc/wp3/judsonb/chinook_hdGBS
conda activate ../hdgbs_env


populations -P "$STACKS_FOLDER" \
        -M "$INFO_FILES_FOLDER"/"$POP_MAP" \
        -t "$NUM_CPU" -p 2 -r 0.6 \
        --ordered-export \
        --fasta-loci \
	--plink \
	--write-random-snp \
        --vcf \
        -O 08_populations
