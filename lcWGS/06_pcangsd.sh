#!/bin/bash

#SBATCH --job-name=lcwgs_pca
#SBATCH --partition=standard
#SBATCH --cpus-per-task=64
#SBATCH --mem-per-cpu=5G
#SBATCH --account=grdi_genarcc
#SBATCH --time=96:00:00

source ~/miniconda3/etc/profile.d/conda.sh
cd /gpfs/fs7/grdi/genarcc/wp3/judsonb/chinook_lcWGS/07e_angsd_m15_maf005_nooutliers # Switch as needed.
conda activate ../../lcwgs_env

export PATH="$HOME/gpfs/fs7/grdi/genarcc/wp3/judsonb/software/pcangsd/bin/:$PATH"

HEADFILE=`ls -1 *.beagle.gz | tail -n 1`
zcat "$HEADFILE" | head -n 1 > header_beagle.txt
cat header_beagle.txt <(gunzip -c *beagle.gz | grep -v "marker") > angsd_c.beagle
gzip -c angsd_c.beagle > angsd_c.beagle.gz

pcangsd -b angsd_c.beagle.gz -t 64 -o ../10_pcangsd/pca_m15_maf005 --maf 0.05





