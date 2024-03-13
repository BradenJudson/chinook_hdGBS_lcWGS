#!/bin/bash

cd /gpfs/fs7/grdi/genarcc/wp3/judsonb/chinook_hdGBS

for i in $(ls -1 02_reads/trimmed_platecat/*.fastq.gz)
do
    basename $i | grep -v "_R2\.fastq\.gz" | perl -pe 's/\.fastq\.gz//'
done > 01_info_files/lane_info.txt

