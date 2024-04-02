#!/bin/bash

#SBATCH --job-name=merge_plates
#SBATCH --partition=standard
#SBATCH --mem-per-cpu=4G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=Braden.Judson@dfo-mpo.gc.ca
#SBATCH --time=01:00:00
#SBATCH --account=grdi_genarcc


cd /gpfs/fs7/grdi/genarcc/wp3/judsonb/chinook_hdGBS/02_reads

# Plate 1, Forward:
cat NS.LH00487_0019.001.D701---B501.Wellband-P1_R1.fastq.gz \
    NS.LH00487_0019.002.D701---B501.Wellband-P1_R1.fastq.gz \
    NS.LH00487_0019.003.D701---B501.Wellband-P1_R1.fastq.gz \
    NS.LH00487_0019.004.D701---B501.Wellband-P1_R1.fastq.gz > trimmed_platecat/trimmed_hdgbs_P1_R1.fastq.gz


# Plate 2, Reverse:
cat NS.LH00487_0019.001.D701---B501.Wellband-P1_R2.fastq.gz \
    NS.LH00487_0019.002.D701---B501.Wellband-P1_R2.fastq.gz \
    NS.LH00487_0019.003.D701---B501.Wellband-P1_R2.fastq.gz \
    NS.LH00487_0019.004.D701---B501.Wellband-P1_R2.fastq.gz > trimmed_platecat/trimmed_hdgbs_P1_R2.fastq.gz


# Plate 2, Forward:
cat NS.LH00487_0019.001.D701---B502.Wellband-P2_R1.fastq.gz \
    NS.LH00487_0019.002.D701---B502.Wellband-P2_R1.fastq.gz \
    NS.LH00487_0019.003.D701---B502.Wellband-P2_R1.fastq.gz \
    NS.LH00487_0019.004.D701---B502.Wellband-P2_R1.fastq.gz > trimmed_platecat/trimmed_hdgbs_P2_R1.fastq.gz

# Plate 2, Reverse:
cat NS.LH00487_0019.001.D701---B502.Wellband-P2_R2.fastq.gz \
    NS.LH00487_0019.002.D701---B502.Wellband-P2_R2.fastq.gz \
    NS.LH00487_0019.003.D701---B502.Wellband-P2_R2.fastq.gz \
    NS.LH00487_0019.004.D701---B502.Wellband-P2_R2.fastq.gz > trimmed_platecat/trimmed_hdgbs_P2_R2.fastq.gz

# Plate 3, Forward:
cat NS.LH00487_0019.001.D702---B501.Wellband-P3_R1.fastq.gz \
    NS.LH00487_0019.002.D702---B501.Wellband-P3_R1.fastq.gz \
    NS.LH00487_0019.003.D702---B501.Wellband-P3_R1.fastq.gz \
    NS.LH00487_0019.004.D702---B501.Wellband-P3_R1.fastq.gz > trimmed_platecat/trimmed_hdgbs_P3_R1.fastq.gz

# Plate 3, Reverse:
cat NS.LH00487_0019.001.D702---B501.Wellband-P3_R2.fastq.gz \
    NS.LH00487_0019.002.D702---B501.Wellband-P3_R2.fastq.gz \
    NS.LH00487_0019.003.D702---B501.Wellband-P3_R2.fastq.gz \
    NS.LH00487_0019.004.D702---B501.Wellband-P3_R2.fastq.gz > trimmed_platecat/trimmed_hdgbs_P3_R2.fastq.gz

# Plate 4, Forward:
cat NS.LH00487_0019.001.D702---B502.Wellband-P4_R1.fastq.gz \
    NS.LH00487_0019.002.D702---B502.Wellband-P4_R1.fastq.gz \
    NS.LH00487_0019.003.D702---B502.Wellband-P4_R1.fastq.gz \
    NS.LH00487_0019.004.D702---B502.Wellband-P4_R1.fastq.gz > trimmed_platecat/trimmed_hdgbs_P4_R1.fastq.gz


# Plate 4, Reverse:
cat NS.LH00487_0019.001.D702---B502.Wellband-P4_R2.fastq.gz \
    NS.LH00487_0019.002.D702---B502.Wellband-P4_R2.fastq.gz \
    NS.LH00487_0019.003.D702---B502.Wellband-P4_R2.fastq.gz \
    NS.LH00487_0019.004.D702---B502.Wellband-P4_R2.fastq.gz > trimmed_platecat/trimmed_hdgbs_P4_R2.fastq.gz

# Plate 5, Forward:
cat NS.LH00487_0019.001.D703---B501.Wellband-P5_R1.fastq.gz \
    NS.LH00487_0019.002.D703---B501.Wellband-P5_R1.fastq.gz \
    NS.LH00487_0019.003.D703---B501.Wellband-P5_R1.fastq.gz \
    NS.LH00487_0019.004.D703---B501.Wellband-P5_R1.fastq.gz > trimmed_platecat/trimmed_hdgbs_P5_R1.fastq.gz

# Plate 5, Reverse:
cat NS.LH00487_0019.001.D703---B501.Wellband-P5_R2.fastq.gz \
    NS.LH00487_0019.002.D703---B501.Wellband-P5_R2.fastq.gz \
    NS.LH00487_0019.003.D703---B501.Wellband-P5_R2.fastq.gz \
    NS.LH00487_0019.004.D703---B501.Wellband-P5_R2.fastq.gz > trimmed_platecat/trimmed_hdgbs_P5_R2.fastq.gz

