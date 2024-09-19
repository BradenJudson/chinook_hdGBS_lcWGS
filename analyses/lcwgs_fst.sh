#!/bin/bash

#SBATCH --account=grdi_genarcc
#SBATCH --job-name=lcwgs_fst
#SBATCH --time=168:00:00
#SBATCH --mem-per-cpu=5G
#SBATCH --cpus-per-task=64
#SBATCH --array=10-1596

# Array is 1/2*(number of populations * (number of populations - 1)).

source ~/miniconda3/etc/profile.d/conda.sh
cd /gpfs/fs7/grdi/genarcc/wp3/judsonb/chinook_lcWGS
conda activate ../lcwgs_env

# Define variables and paths.
GENOMEFOLDER=/gpfs/fs7/grdi/genarcc/common/genomes/Chinook_Salmon/index_faidx
GENOME="GCF_018296145.1_Otsh_v2.0_genomic.fna"
POPLIST=`ls -1 13_pop_bams/text_files | sed -e 's/_bams.txt//g'`
OUTFOLDER="14_fst/sfs_outputs"
SITES="01_info_files/chinook_maf005_8MSNPs.txt"
TDSFSDIR="14_fst/2dsfs"
FSTINDEX="14_fst/fst_index"
GFST="14_fst/global_fst"
NCORES=64

# For each population, estimate the site frequency spectrum (SFS) using -doSaf 1.
for POP in $POPLIST
do
        echo $POP
        angsd -b 13_pop_bams/text_files/"$POP"_bams.txt \
        -ref "$GENOMEFOLDER"/"$GENOME" \
        -anc "$GENOMEFOLDER"/"$GENOME" \
        -out "$OUTFOLDER"/"$POP" \
        -GL 1 -doSaf 1 -doCounts 1 \
        -sites "$SITES"

done

# Make a list of all population names and a separate file of all possible pairs of populations.
ls -1 13_pop_bams/text_files | sed -e 's/_bams.txt//g' > 01_info_files/pops.txt
awk '{T[NR]=$1} END {for (i=1; i<NR; i++) {for (j=i+1; j<=NR; j++) print T[i], T[j]}}' 01_info_files/pops.txt > 01_info_files/population_pairs.txt

# Isolate one pair of populations at a time.
POP1=`cat 01_info_files/population_pairs.txt | head -n $SLURM_ARRAY_TASK_ID | tail -n 1 | cut -d ' ' -f 1`
POP2=`cat 01_info_files/population_pairs.txt | head -n $SLURM_ARRAY_TASK_ID | tail -n 1 | cut -d ' ' -f 2`

# Calculate the (folded) 2dSFS prior.
realSFS "$OUTFOLDER"/"$POP1".saf.idx "$OUTFOLDER"/"$POP2".saf.idx -fold 1 > "$TDSFSDIR"/"$POP1"_"$POP2".sfs

# Prepare FST file. Make sure to specify the folded SFS file from earlier.
realSFS fst index "$OUTFOLDER"/"$POP1".saf.idx "$OUTFOLDER"/"$POP2".saf.idx -sfs "$TDSFSDIR"/"$POP1"_"$POP2".sfs -fstout "$FSTINDEX"/"$POP1"_"$POP2" -P "$NCORES"

# Global FST estimates.
realSFS fst stats "$FSTINDEX"/"$POP1"_"$POP2".fst.idx -P "$NCORES" > "$GFST"/"$POP1"_"$POP2".globalFst

# Merge all pairwise global Fst estimates into a single longform file. Subsequent manipulation into a pairwise matrix done in R.
awk '{print FILENAME "\t", $0}' "$GFST"/*.globalFst | sed '1i file\tFstUnweighted\tFstWeighted' - > "$GFST"/all_fst_estimates.txt