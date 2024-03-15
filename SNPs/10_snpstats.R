setwd("~/ots_landscape_genetics/SNPs")

library(tidyverse); library(ggplot2); library(karyoploteR)

# Read in genome organization stats (from NCBI).
genome <- read.table("./info_files/ch_genome_seqrep.txt", 
                     sep = "\t", header = T) %>%
  # Retain autosomes only for present purposes.
  filter(Role != "unplaced-scaffold" & Molecule.type == "Chromosome") %>% 
  # Reorganize for reading into toGRanges (see ?toGRanges).
  select(c("Chromosome.name", "Seq.length")) %>% 
  `colnames<-`(., c("chr", "end")) %>% 
  mutate(chr = paste0("Ots", str_pad(chr, 
               width=2, side = "left", pad = "0")),
         start = 1) %>% 
  relocate(start, .after = chr)

snps <- toGRanges() ### ADD SNPS HERE

tiff(filename = "./stats/chrom_snps.tiff", width = 10, 
     height = 15, res = 300, units = "in")

# Save the tiff file here.
(p <- plotKaryotype(genome = toGRanges(genome))); dev.off()

