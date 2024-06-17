#!/bin/bash

#SBATCH --job-name=lcwgs_clipov
#SBATCH --partition=standard
#SBATCH --mem-per-cpu=5G
#SBATCH --ntasks=1
#SBATCH --mail-user=Braden.Judson@dfo-mpo.gc.ca
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=64
#SBATCH --comment="registry.maze.science.gc.ca/ssc-hpcs/generic-job:ubuntu22.04"
#SBATCH --account=grdi_genarcc


# Source conda environment.
source ~/.bashrc

INPUT="05_deduplicated_bams"
OUTPUT="06_clipped_bams"

# Remove alignment overlaps within paired reads

for file in $(ls "$INPUT"/*.bam | perl -pe 's/\.bam//g')
do
        name=$(basename "$file")
        echo "Clipping sample: $name"

        bam clipOverlap \
                --in "$INPUT"/"$name".bam \
                --out "$OUTPUT"/"$name".clip.bam \
                --stats

        samtools index "$OUTPUT"/"$name".clip.bam

done