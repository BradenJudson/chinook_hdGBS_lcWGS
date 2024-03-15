#!/bin/bash

#SBATCH --job-name=hdgbs_gstacks
#SBATCH --partition=standard
#SBATCH --mail-user=Braden.Judson@dfo-mpo.gc.ca
#SBATCH --time=60:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5G
#SBATCH --ntasks=1
#SBATCH --account=grdi_genarcc

# Source conda environment.
source ~/miniconda3/etc/profile.d/conda.sh
cd /gpfs/fs7/grdi/genarcc/wp3/judsonb/chinook_hdGBS
conda activate ../hdgbs_env

STACKS_FOLDER="06_stacks"
SAMPLE_FOLDER="05_alignment"
INFO_FILES="01_info_files"
POP_MAP="population_map_filt.txt"

gstacks -I "$SAMPLE_FOLDER" -S ".1.sorted.bam" \
        -M "$INFO_FILES"/"$POP_MAP" \
        -O "$STACKS_FOLDER" \
        --max-clipped 0.1 \
        --rm-pcr-duplicates