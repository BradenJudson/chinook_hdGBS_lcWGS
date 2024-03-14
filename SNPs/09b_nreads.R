library(tidyverse); library(ggplot2)
library(ggExtra); library(cowplot)

setwd("~/ots_landscape_genetics/SNPs")

Th <- 1.5e6 # Reads threshold for pass/fail.

# This file contains the read numbers per bam file.
reads <- read.table(file = "stats/sample_reads.txt",
                  col.names = c("Sample", "Reads")) %>% 
  filter(Reads > 0) %>% 
  mutate(Sample = gsub('.1.sorted.bam', "", Sample),
         pass = case_when(Reads >= Th ~ "Y",
                          Reads <  Th ~ "N"))

write.csv(reads, "stats/sample_reads.csv")


# Merge population-level info with read data.
sample_info <- read.csv("info_files/hdGBS_sampleinfo.csv") %>% 
  merge(., reads, by = "Sample")

(RP <- ggplot(data = sample_info, 
       aes(x = Population, y = Reads/1e6)) + 
  geom_hline(yintercept = Th/1e6, colour = "blue2",
             linetype ="dashed") +
  geom_hline(yintercept = mean(sample_info$Reads/1e6),
             color = "blue", linetype = "dashed") +
  geom_boxplot(alpha = 1/10) +
  geom_point(size  = 3/2, 
             shape = 21,
             color = "black",
             aes(fill = pass)) +
  scale_fill_manual(values = c("red1", "gray80")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90,
                      vjust = 0.5, hjust = 1),
        legend.position = "none") +
  labs(x = NULL, y = "Reads (millions)"))

(RPH <- ggMarginal(RP, type = "histogram", margins = "y",
           binwidth = 1, fill = "gray90", size = 10))

# Use cowplot here to save multi-plot configuration more easily.
save_plot("./stats/pop_reads.tiff", RPH, ncol = 2, base_height = 6)

# Summary stats overall.
summary(sample_info$Reads)
quantile(sample_info$Reads, probs = c(0.01, 0.05, 0.08, 0.1, 0.9, 0.95, 0.99))

# Summary stats by population.
(pop_reads <- sample_info %>% 
  filter(pass == "Y") %>% 
  group_by(Population) %>% 
  summarise(mean = mean(Reads, na.rm = T),
            sD   = sd(Reads,   na.rm = T),
            min  = min(Reads,  na.rm = T),
            max  = max(Reads,  na.rm = T),
            n    = n()))

hist(pop_reads$n); summary(pop_reads$n)

write.csv(pop_reads, "stats/pop_reads.csv", row.names = F)

# Identify fail-rates per population.
(pop_fails <- sample_info[sample_info$pass == "N",] %>% 
  group_by(Population) %>% 
  tally())

sum(pop_fails$n)
sum(pop_fails$n)/nrow(sample_info)*100

write.csv(pop_fails, "stats/pop_fails.csv", row.names = F)


