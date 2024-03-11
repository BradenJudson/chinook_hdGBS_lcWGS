#!/bin/bash

cd /gpfs/fs7/grdi/genarcc/wp3/judsonb/chinook_hdGBS

# Following section 4.3.2 here:
# https://catchenlab.life.illinois.edu/stacks/manual/#popmap

grep -v "#" 01_info_files/hdGBS_sampleinfo.csv | \
        cut -d ',' -f 4,5 | \
        sed 's/\s.*$//' | \
        tr ',' '\t' | \
        sort -u > 01_info_files/population_map.txt
