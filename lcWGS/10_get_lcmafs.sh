#!/bin/bash

cd /gpfs/fs7/grdi/genarcc/wp3/judsonb/chinook_lcWGS/07e_angsd_m15_maf005_nooutliers

# Pick first mafs file and pull header row.
zcat chinook_angsd_NC_056450.1:10000000-20000000.mafs.gz | head -n 1 > header.txt

# Add header above to the concatenation of every mafs file in the directory.
# Tail command removes their headers and awk subsets relevant columns.
cat header.txt <(cat *.mafs.gz | gunzip -c | tail -n +2 ) |\
        awk '{print $1,$2,$6}' - > ../15_lc_afs/lcafs_maf005.mafs