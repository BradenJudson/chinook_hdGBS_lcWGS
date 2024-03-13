#!/bin/bash


source ~/miniconda3/etc/profile.d/conda.sh
cd /gpfs/fs7/grdi/genarcc/wp3/judsonb/chinook_hdGBS/05_alignment
conda activate ../../hdgbs_env


for file in *.sorted.bam;
do
        sample="$file"
        nreads=$(samtools view -c -F 260 "$file");
        echo ${sample}  ${nreads}
done > ../sample_reads.txt
