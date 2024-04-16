#!/bin/bash


#SBATCH --job-name=hdgbs_bwa
#SBATCH --partition=standard
#SBATCH --mem-per-cpu=5G
#SBATCH --ntasks=1
#SBATCH --mail-user=Braden.Judson@dfo-mpo.gc.ca
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=64
#SBATCH --account=grdi_genarcc


# Global variables
GENOMEFOLDER=/gpfs/fs7/grdi/genarcc/common/genomes/Chinook_Salmon
GENOME="GCF_018296145.1_Otsh_v2.0_genomic.fna"
DATAFOLDER="04_all_samples"
NCPU=1
ALIGNFOLDER="05_alignment"


# Source conda environment.
source ~/.profile


for file in $(ls -1 "$DATAFOLDER"/*.1.fq.gz)
do

        file2=$(echo "$file" | perl -pe 's/\.1.fq.gz/\.2.fq.gz/')
        echo "Aligning $file $file2"

        name=$(basename "$file")
        name2=$(basename "$file2")

        bwa mem -t "$NCPU" \
                "$GENOMEFOLDER"/"$GENOME" "$DATAFOLDER"/"$name" "$DATAFOLDER"/"$name2" 2> /dev/null |
                samtools view -Sb -q 20 - > "$ALIGNFOLDER"/"${name%.fq.gz}".bam

        samtools sort --threads "$NCPU" -o "$ALIGNFOLDER"/"${name%.fq.gz}".sorted.bam \
                "$ALIGNFOLDER"/"${name%.fq.gz}".bam

        samtools index "$ALIGNFOLDER"/"${name%.fq.gz}".sorted.bam

done
