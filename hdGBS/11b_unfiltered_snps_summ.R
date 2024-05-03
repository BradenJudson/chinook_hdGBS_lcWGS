setwd("~/ots_landscape_genetics/hdGBS")

library(tidyverse); library(ggplot2); library(ggExtra); library(cowplot)


################################################################################
# SNP Depths -------------------------------------------------------------------
################################################################################

SNPdepth <- read.table("./stats/unfiltered_stats_complete/out.ldepth.mean", header = T); colnames(SNPdepth)

(mean_snp_depth <- summary(SNPdepth$MEAN_DEPTH)); sd(SNPdepth$MEAN_DEPTH)
quantile(SNPdepth$MEAN_DEPTH, probs = c(0.1, 0.5, 0.9, 0.95, 0.99, 0.999, 0.9999))
nrow(SNPdepth[SNPdepth$MEAN_DEPTH > 100, ])

# (rsd <- sd(SNPdepth$MEAN_DEPTH))
# (dep_thr <- mean(SNPdepth$MEAN_DEPTH)+c(rsd*-1, rsd, rsd*2))

# DEPTHLT <- 3
# DEPTHUT <- 13
# 
# SNPdepth$filter <- case_when(SNPdepth$MEAN_DEPTH <  DEPTHLT | SNPdepth$MEAN_DEPTH >  DEPTHUT ~ "F",
#                              SNPdepth$MEAN_DEPTH >= DEPTHLT | SNPdepth$MEAN_DEPTH <= DEPTHUT ~ "P")

(snp_depthhist <- ggplot() +
  geom_histogram(data = SNPdepth,
                 aes(MEAN_DEPTH), bins = 150, 
                 color = "black", fill = "gray80") +
  theme_bw() + xlim(3, 30) +
  geom_vline(xintercept = DEPTHLT, color = "red") +
  labs(x = "Mean SNP depth", y = "Frequency"))

ggsave("./stats/snp_depth.tiff", width = 12,
       height = 6, dpi = 300)

SNP_d3 <- SNPdepth[SNPdepth$MEAN_DEPTH >= 3 & SNPdepth$MEAN_DEPTH < 100, ]

# SNP Missingness --------------------------------------------------------------

SNPmiss <- read.table("./stats/unfiltered_stats_complete/out.lmiss", header = T)


MISST <- 4/10
SNPmiss$filter <- case_when(SNPmiss$F_MISS >= MISST ~ "F",
                            SNPmiss$F_MISS <  MISST ~ "P")

(snp_miss_summ <- summary(SNPmiss$F_MISS))

(snp_misshist <- ggplot() +
    geom_histogram(data = SNPmiss, 
                   aes(F_MISS), bins = 180, color = "black", 
                   fill = "gray80")  +
    labs(x = "SNP Missingness (%)", y = "Frequency") +
    theme_bw())


ggsave("./stats/snp_missingnessh.tiff", width = 12,
       height = 6, dpi = 300)

SNP_lowmiss <- SNPmiss[SNPmiss$F_MISS < 0.4, ]; nrow(SNP_lowmiss)

# MAF --------------------------------------------------------------------------

SNPhet <- read.table("./stats/unfiltered_stats_complete/out.frq", header = T, row.names = NULL) %>% 
  `colnames<-`(., c("CHROM", "POS", "NALLELES", "N_CHR", "MajAll", "MinAll")) %>% 
  mutate(MAF = as.numeric(sub(".*:", "", MinAll)))

# Minor allele frequency cutoff of 2%
MAFT <- 0.02

SNPhet$filter <- case_when(SNPhet$MAF >= MAFT ~ "P",
                           SNPhet$MAF <  MAFT ~ "F")

(maf_summ <- summary(SNPhet$MAF))

(mafhist <- ggplot() +
    geom_histogram(data = SNPhet, aes(MAF),
                   bins = 250, color = "black", fill = "gray80") +
    labs(x = "Minor Allele Frequency", y = "Frequency") +
    geom_vline(xintercept = MAFT, color = "red") +
    scale_x_continuous(limits = c(0, 0.5)) +
    theme_bw())

ggsave("./stats/snp_MAF.tiff", width = 12,
       height = 6, dpi = 300)

maf5 <- SNPhet[SNPhet$MAF > MAFT, ]; nrow(maf5)
write.table(maf5[,c("CHROM", "POS")], "./filters/maf5.tsv", 
            quote = F, row.names = F, col.names = F)



SNP_count <- read.table("./stats/unfiltered_stats_complete/out.frq.count", header = T, row.names = NULL) %>% 
  mutate(MA = as.numeric(sub(".*:", "", X.ALLELE.COUNT.))) %>% 
  filter(MA > 0)


(mac <- ggplot(data = SNP_count, aes(MA)) +
  geom_histogram(bins = 250, color = "black", fill = "gray80") +
    theme_bw() + labs(x = "Minor allele count", y = "Frequency") +
    xlim(0, 50)) 

ggsave("./stats/mac.tiff", width = 12,
       height = 6, dpi = 300)

# Write Snp stats --------------------------------------------------------------
(snp_stats <- rbind(maf_summ, snp_miss_summ, mean_snp_depth))
write.csv(snp_stats, "./stats/snps_stats.csv")

(qplots <- cowplot::plot_grid(mafhist, snp_misshist, snp_depthhist, 
                              indv_readhist, ncol = 2))

ggsave("./stats/genotype_filtering_2x2.tiff", width = 12, height = 8, dpi = 300)


################################################################################
# Individual Read Depth --------------------------------------------------------
################################################################################


# Read in VCFtools output for missing data per individual. 
depth <- read.table("./stats/unfiltered_stats_complete/out.idepth", header = T)

# Plot missing data histogram.
ggplot(data = depth, aes(MEAN_DEPTH)) +
  geom_histogram(bins = 150,
                 color = "black",
                 fill = "gray80") +
  labs(x = "Mean depth",
       y = "Frequency") +
  theme_bw() +
  scale_x_continuous(breaks = seq(0, 40, 2))

ggsave("./stats/indv_depth.tiff", width = 12,
       height = 6, dpi = 300)

# Summarize individual missing data and select threshold for discarding data.
(mean_depth <- summary(depth$MEAN_DEPTH))
depth.threshold <- 300e3
sum(depth$N_SITES < depth.threshold)

# Select individuals passing depth filtering step.
depth.keep <- depth[depth$N_SITES > depth.threshold, "INDV"] 


# INDV Missingness -------------------------------------------------------------

imiss <- read.table("./stats/unfiltered_stats_complete/out2.imiss", header = T)
(indv_miss_p <- summary(imiss$F_MISS))
(indv_miss_n <- summary(imiss$N_MISS))

quantile(imiss$F_MISS, probs = c(0.8, 0.9, 0.95, 0.99))

nrow(imiss[imiss$F_MISS < 0.90,])

(imisshist <- ggplot(data = imiss, aes(F_MISS)) +
    geom_histogram(bins = 150, color = "black",
                   fill = "gray80") +
    labs(y = "Frequency", x = "Missing positions (% Total)") +
    theme_bw())

ggsave("./stats/indv_missing_n.tiff", width = 12,
       height = 6, dpi = 300)

indmiss <- imiss[imiss$F_MISS > 0.3, ] 


# INDV reads -------------------------------------------------------------------

ThU <- 2e6 # Reads threshold for pass/fail [[arbitary]]
ThL <- 130e6

# This file contains the read numbers per bam file.
reads <- read.table(file = "stats/sample_reads2.txt",
                    col.names = c("Sample", "Reads")) %>% 
  filter(Reads > 0) %>% 
  mutate(Sample = gsub('.1.sorted.bam', "", Sample),
         pass = case_when(Reads >= ThU & Reads < ThL  ~ "Y",
                          Reads <  Th | Reads  > ThU  ~ "N" ))

quantile(reads$Reads, probs = c(.01, .05, .10, .5, .9, .95, .99))
mean(reads$Reads); sd(reads$Reads)

# Merge population-level info with read data.
sample_info <- read.csv("info_files/hdGBS_sampleinfo.csv") %>% 
  merge(., reads, by = "Sample")

# Summary stats by population.
(pop_reads <- sample_info %>%  
    group_by(Population) %>% 
    summarise(mean = mean(Reads, na.rm = T),
              sD   = sd(Reads,   na.rm = T),
              min  = min(Reads,  na.rm = T),
              max  = max(Reads,  na.rm = T),
              n    = n()))

write.csv(pop_reads, "stats/pop_reads.csv", row.names = F)

(RP <- ggplot(data = sample_info, 
              aes(x = Population, y = Reads/1e6)) + 
    geom_hline(yintercept = mean(sample_info$Reads/1e6)-sd(sample_info$Reads/1e6), 
               colour = "blue2", linetype ="dashed", alpha = 1/4) +
    geom_hline(yintercept = mean(sample_info$Reads/1e6)+sd(sample_info$Reads/1e6), 
               colour = "blue2", linetype ="dashed", alpha = 1/4) +
    geom_hline(yintercept = mean(sample_info$Reads/1e6),
               color = "blue", linetype = "dashed") +
    geom_boxplot(alpha = 1/10) +
    geom_point(size  = 3/2, 
               shape = 21,
               color = "black",
               fill = 'grey80') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90,
                                     vjust = 0.5, hjust = 1),
          legend.position = "none") +
    labs(x = NULL, y = "Reads (millions)") +
    scale_y_continuous(breaks = c(0, 5, 25, 50, 100, 150, 200),
                       expand = c(0.01, 0.01)))

# Puts marginal histogram on right-most y-axis.
# Shows skewed distribution.
(RPH <- ggMarginal(RP, type = "histogram", margins = "y",
                   binwidth = 2, fill = "gray90", size = 10))

# Use cowplot here to save multi-plot configuration more easily.
save_plot("./stats/pop_reads_1.tiff", RPH, ncol = 2, base_height = 6)

(indv_readhist <- ggplot() +
  geom_histogram(data = sample_info[sample_info$pass == "Y", ], 
                 aes(Reads/1e6), color = "black", fill = "gray80", bins = 150) +
  geom_histogram(data = sample_info[sample_info$pass == "N", ], 
                 aes(Reads/1e6), color = "gray80", fill = "gray99", bins = 150) +
  theme_bw() +  geom_vline(xintercept = Th/1e6, color = "red") +
  labs(x = "Reads per sample (millions)", y = "Frequency") +
  scale_x_continuous(breaks = seq(0, 200, 20)) +
  scale_y_continuous(breaks = seq(0, 30, 5)))

ggsave("./stats/read_depth_hist.tiff", width = 12,
       height = 6, dpi = 300)

ind_fails <- sample_info[sample_info$pass == "N", "Sample"]

write.table(ind_fails, "stats/ind_fails.txt", quote = FALSE, 
            row.names = FALSE, col.names = FALSE, sep = "\t")

