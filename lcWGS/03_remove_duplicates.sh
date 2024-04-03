#!/bin/bash

#SBATCH --job-name=lcwgs_dedup
#SBATCH --partition=standard
#SBATCH --mem-per-cpu=5G
#SBATCH --mail-user=Braden.Judson@dfo-mpo.gc.ca
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=32
#SBATCH --account=grdi_genarcc


source ~/miniconda3/etc/profile.d/conda.sh
cd /gpfs/fs7/grdi/genarcc/wp3/judsonb/chinook_lcWGS
conda activate ../lcwgs_env


INPUT="04_all_alignments"
OUTPUT="05_deduplicated_bams"
METRICS="11_metrics"
JAVA_OPTS="-Xmx80G"
TMPDIR="99_tmp"


# Remove duplicate alignments
for file in $(ls "$INPUT"/*.bam | perl -pe 's/\.bam//g')
do
        name=$(basename "$file")
        echo "Deduplicating sample: $name"

        picard $JAVA_OPTS MarkDuplicates \
                I="$INPUT"/"$name".bam \
                O="$OUTPUT"/"$name".dedup.bam \
                M="$METRICS"/"$name".metrics.txt \
                TMP_DIR="$TMPDIR" \
                REMOVE_DUPLICATES=true

done
