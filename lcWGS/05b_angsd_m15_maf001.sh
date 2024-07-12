#!/bin/bash

#SBATCH --job-name=ch_angsd
#SBATCH --account=grdi_genarcc
#SBATCH --partition=standard
#SBATCH --output=logs_may23/angsd_chr.%A_%a.out
#SBATCH --error=logs_may23/angsd_chr.%A_%a.err
#SBATCH --time=14:00:00
#SBATCH --nodes=1
#SBATCH --mem=110GB
#SBATCH --cpus-per-task=1
#SBATCH --array=1-217

source ~/miniconda3/etc/profile.d/conda.sh
cd /gpfs/fs7/grdi/genarcc/wp3/judsonb/chinook_lcWGS
conda activate ../lcwgs_env

REGION=`cat ../chinook_10mb_genome.txt | head -n $SLURM_ARRAY_TASK_ID | tail -n 1`
GENOMEFOLDER=/gpfs/fs7/grdi/genarcc/common/genomes/Chinook_Salmon/index_faidx
GENOME="GCF_018296145.1_Otsh_v2.0_genomic.fna"
BAMFILES="bam_list_n453.txt"
OUTFOLDER="07f_angsd_maf_maf001_nooutliers"

#ls -1 06_clipped_bams/*.bam > 01_info_files/"$BAMFILES"

angsd -b 01_info_files/"$BAMFILES" -ref "$GENOMEFOLDER"/"$GENOME" \
        -out "$OUTFOLDER"/chinook_angsd_"$REGION" \
        -nThreads 8 \
        -r "$REGION" \
        -uniqueOnly 1 \
        -remove_bads 1 \
        -SNP_pval 1e-10  \
        -minMapQ 30 \
        -minMaf 0.01 \
        -minQ 30 \
        -setMinDepth 453 \
        -setMaxDepth 2000 \
        -minInd 385 \
        -doMaf 1 \
        -doMajorMinor 1 \
        -doCounts 1 \
        -GL 1 -doGLF 2 \
        -only_proper_pairs 1 \
        -dumpCounts 2 \
        -rmTriallelic 1e-6 \
        -doBcf 1 -doPost 1 \
        -doGeno 5 -doQsDist 1