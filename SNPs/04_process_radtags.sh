#!/bin/bash


#SBATCH --job-name=hdgbs_pradtags
#SBATCH --partition=standard
#SBATCH --mem-per-cpu=5G
#SBATCH --ntasks=1
#SBATCH --mail-user=Braden.Judson@dfo-mpo.gc.ca
#SBATCH --time=16:00:00
#SBATCH --cpus-per-task=64
#SBATCH --account=grdi_genarcc


source ~/miniconda3/etc/profile.d/conda.sh
cd /gpfs/fs7/grdi/genarcc/wp3/judsonb/chinook_hdGBS
conda activate ../hdgbs_env


ENZYME="bfaI"
TRIM_LENGTH=120
NCPU=64


# Load GNU parallel.
. ssmuse-sh -x main/opt/parallel/parallel-202109022

export stacks
export process_radtags

INFO_FILES="01_info_files"

# Pipe works such that the lane info strings are read as the third argument in the utility script.
# I.e., in process_radtags_1_enzyme_pe.sh, $3 will denote the lane being processed.
cat $INFO_FILES/lane_info.txt |
        parallel -j $NCPU ./00_scripts/utility_scripts/process_radtags_1_enzyme_pe.sh "$TRIM_LENGTH" "$ENZYME"
