#!/bin/bash


#SBATCH --job-name=lcwgs_bwa3
#SBATCH --partition=standard
#SBATCH --mem-per-cpu=5G
#SBATCH --ntasks=1
#SBATCH --mail-user=Braden.Judson@dfo-mpo.gc.ca
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=64
#SBATCH --account=grdi_genarcc


# Source conda environment.
source ~/miniconda3/etc/profile.d/conda.sh
cd /gpfs/fs7/grdi/genarcc/wp3/judsonb/chinook_lcWGS
conda activate ../lcwgs_env


GENOME_FOLDER=/gpfs/fs7/grdi/genarcc/common/genomes/Chinook_Salmon
GENOME="GCF_018296145.1_Otsh_v2.0_genomic.fna"
TRIMMED_FOLDER="03d_trimmed_reads_alignment4"
ALIGNED_FOLDER="04_all_alignments"
NCPU=64
TEMP_FOLDER="99_tmp"
SAMPLE_FILE="samples.txt"


cat "$SAMPLE_FILE" |
while read file
do
        base=$(basename $file | perl -pe 's/_R1.fastq.gz//')
        echo "Aligning $base"

        ID=$(echo $base | perl -pe 's/\.[a-zA-Z]*-.*$//g')
        PU=$(echo $base | perl -pe 's/NS\.LH00487_//g' | perl -pe 's/\.IDT_i7.*$//g')
        LB=$(echo $base | perl -pe 's/^.*[0-9]*\.IDT/IDT/g')
        SM=$(echo $base | perl -pe 's/NS.*i5.*_[0-9]*\.//g')

        bwa mem -t "$NCPU" \
               "$GENOME_FOLDER"/"$GENOME" \
                "$TRIMMED_FOLDER"/"$base"_R1.fastq.gz \
                "$TRIMMED_FOLDER"/"$base"_R2.fastq.gz \
                -R "@RG\tID:$ID\tPL:illumina\tPU:$PU\tSM:$SM\tLB:$LB" |
                samtools view -b -q 20 - |
                samtools sort -T "$TEMP_FOLDER"/"$base" - > "$ALIGNED_FOLDER"/"$base".bam

        samtools index "$ALIGNED_FOLDER"/"$base".bam

        samtools flagstat "$ALIGNED_FOLDER"/"$base".bam > "$ALIGNED_FOLDER"/"$base".flagstat

        samtools idxstat "$ALIGNED_FOLDER"/"$base".bam > "$ALIGNED_FOLDER"/"$base".idxstat
done
