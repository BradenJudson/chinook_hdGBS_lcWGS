setwd("~/ots_landscape_genetics/hdGBS")

library(tidyverse); library(ggplot2); library(karyoploteR)

# Read in genome organization stats (from NCBI).
genome <- read.table("./info_files/ch_genome_seqrep.txt", 
                     sep = "\t", header = T) %>%
  # Retain autosomes only for present purposes.
  filter(Role != "unplaced-scaffold" & Molecule.type == "Chromosome") %>% 
  # Reorganize for reading into toGRanges (see ?toGRanges).
  select(c("Chromosome.name", "Seq.length")) %>% 
  `colnames<-`(., c("chr", "end")) %>% 
  # Relabel chromosomes to a consistent "Ots##" format.
  mutate(chr = paste0("Ots", str_pad(chr, 
               width=2, side = "left", pad = "0")),
         start = 1,
         tpos = cumsum(end)) %>% 
  relocate(start, .after = chr)

# Plot base autosomes.
(base <- ggplot() +
    theme_minimal() +
    theme_light() +
       geom_bar(data = genome,
             aes(x = chr, y = end/1e6),
             stat = 'identity',
             fill = 'grey80',
             width = 1/5) +
    labs(x = NULL, y = "Base pair position (Mbp)") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          axis.ticks.x = element_blank(),
          panel.grid = element_blank()))

ggsave("./stats/chrom_snps.tiff", width = 15, height = 10, dpi = 300)


# Read in SNPs.



# Depth -------------------------------------------------------------------

depth <- read.table("./stats/out.idepth.txt", header = T)

ggplot(data = depth, aes(N_SITES/1e3)) +
  geom_histogram(bins = 100,
                 color = "black",
                 fill = "gray80") +
  labs(x = "Sites (Thousands)",
       y = "Frequency") +
  theme_bw() +
  scale_x_continuous(breaks = seq(0, 1600, 200))

ggsave("./stats/depth.tiff", width = 12,
       height = 6, dpi = 300)

summary(depth$N_SITES)
sum(depth$N_SITES < 300e3)
