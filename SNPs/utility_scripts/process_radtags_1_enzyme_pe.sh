#!/bin/bash


# Source conda environment.
source ~/miniconda3/etc/profile.d/conda.sh
cd /gpfs/fs7/grdi/genarcc/wp3/judsonb/chinook_hdGBS
conda activate ../hdgbs_env


TRIM_LENGTH=$1
ENZYME1=$2
LANE=$3
INFO_FILES="01_info_files"
TEMPBARCODES=.temp.${LANE}.barcodes


# Create temporary barcode files specific to each lane/plate.
grep -vE "^$" "$INFO_FILES"/hdGBS_sampleinfo.csv | \
    grep "$LANE" | \
    grep -v "Barcode" | \
    cut -d ',' -f 2 > "$INFO_FILES"/"$TEMPBARCODES"


# Demultiplex samples.
# q = quality, remove low-quality samples.
# b = directory housing barcode files produced above.
# c = clean, remove reads with uncalled bases.
# r = rescue barcodes and RAD-Tag cut sites.
# barcode-dist-1: Number of allowed mismatches when rescuing single-end barcodes.
process_radtags -i gzfastq \
        -1 02_reads/trimmed/"$LANE"".fastq.gz" \
        -2 02_reads/trimmed/$(echo "$LANE" | perl -pe 's/_R1/_R2/')".fastq.gz" \
        -o 03_samples/"$LANE" \
        -b "$INFO_FILES"/"$TEMPBARCODES" \
        -c -r -t "$TRIM_LENGTH" \
        -q \
        -s 0 \
        -E phred33 \
        --barcode-dist-1 2 \
        --renz_1 "$ENZYME1"
