#!/bin/bash


# https://github.com/enormandeau/stacks_workflow/blob/master/00-scripts/03_rename_samples_pe.sh
# ^ Slight modification of this.


INFO_FILES="01_info_files"
ALL_SAMPLES_FOLDER="04_all_samples"
SAMPLES_FOLDER="03_samples"

cd /gpfs/fs7/grdi/genarcc/wp3/judsonb/chinook_hdGBS

# Rename left files.
grep -vE '^#' "$INFO_FILES"/hdGBS_sampleinfo.csv |
    cut -d ',' -f 1,2,4 |
    perl -pe 's/\.f(ast)*q\.gz//' |
    perl -pe 's/,/\/sample_/' |
    perl -pe 's/([ACTG]+),/\1.fq.gz\t/' |
    awk '{print $1"\t04_all_samples/"$2".fq.gz"}' |
    perl -pe 's/\.fq\.gz/.1.fq.gz/g' > renaming_01l.txt


cut -f 2 renaming_01l.txt | sort -u > renaming_02l.txt

# Rename right files.
grep -vE '^#' "$INFO_FILES"/hdGBS_sampleinfo.csv |
    cut -d ',' -f 1,2,4 |
    perl -pe 's/\.f(ast)*q\.gz//' |
    perl -pe 's/,/\/sample_/' |
    perl -pe 's/([ACTG]+),/\1.fq.gz\t/' |
    awk '{print $1"\t04_all_samples/"$2".fq.gz"}' |
    perl -pe 's/\.fq\.gz/.2.fq.gz/g' > renaming_01r.txt

cut -f 2 renaming_01r.txt | sort -u > renaming_02r.txt

cat renaming_02l.txt |
    while read -r i
    do
        echo Treating: "$i"
        rm $ALL_SAMPLES_FOLDER/"$i" 2> /dev/null
        num_copies=$(grep "$i" renaming_01l.txt | cut -f 1 | wc -l)
        echo -n "  Sample found $num_copies times: "
        if [[ $num_copies -eq 1 ]]
        then
            echo "Creating a link"
            grep "$i" renaming_01l.txt |
                cut -f 1 |
                while read -r j
                do
                    echo "    Linking: $SAMPLES_FOLDER/""$j"
                    ln $SAMPLES_FOLDER/"$j" "$i"
                done
        else
            echo "Merging samples"
            grep "$i" renaming_01l.txt |
                cut -f 1 |
                while read -r j
                do 
                    echo "    Copying: $SAMPLES_FOLDER/""$j"
                    cat $SAMPLES_FOLDER/"$j" >> "$i"
                done
        fi
    done

cat renaming_02r.txt |
    while read -r i
    do
        echo Treating: "$i"
        rm $ALL_SAMPLES_FOLDER/"$i" 2> /dev/null
        num_copies=$(grep "$i" renaming_01r.txt | cut -f 1 | wc -l)
        echo -n "  Sample found $num_copies times: "
        if [[ $num_copies -eq 1 ]]
        then
            echo "Creating a link"
            grep "$i" renaming_01r.txt |
                cut -f 1 |
                while read -r j
                do
                    echo "    Linking: $SAMPLES_FOLDER/""$j"
                    ln $SAMPLES_FOLDER/"$j" "$i"
                done
        else
            echo "Merging samples"
            grep "$i" renaming_01r.txt |
                cut -f 1 |
                while read -r j
                do
                    echo "    Copying: $SAMPLES_FOLDER/""$j"
                    cat $SAMPLES_FOLDER/"$j" >> "$i"
                done
        fi
    done
