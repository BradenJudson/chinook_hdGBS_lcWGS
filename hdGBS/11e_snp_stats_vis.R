setwd("~/ots_landscape_genetics/hdGBS/stats/filtered_stats")

library(tidyverse); library(ggplot2)


################################################################################
# indv depth -------------------------------------------------------------------
################################################################################


idepth <- read.table("out.idepth", header = T); colnames(SNPdepth)

summary(idepth$MEAN_DEPTH)

ggplot(data = idepth, aes(MEAN_DEPTH)) +
  geom_histogram(bins = 150,
                 color = "black",
                 fill = "gray80") +
  labs(x = "Mean depth",
       y = "Frequency") +
  theme_bw() +
  scale_x_continuous(breaks = seq(0, 21, 1))

ggsave("indv_depth.tiff", width = 12, 
       height = 6, dpi = 300)


################################################################################
# indv miss --------------------------------------------------------------------
################################################################################

imiss <- read.table("out.imiss", header = T)

(indv_miss_p <- summary(imiss$F_MISS))
(indv_miss_n <- summary(imiss$N_MISS))

(imisshist <- ggplot(data = imiss, aes(F_MISS)) +
    geom_histogram(bins = 150, color = "black",
                   fill = "gray80") +
    labs(y = "Frequency", x = "Missing positions (% Total)") +
    theme_bw())

nrow(imiss[imiss$F_MISS > 0.6, ]); nrow(imiss)

ggsave("indv_missing_n.tiff", width = 12, 
       height = 6, dpi = 300)

write.table(imiss[imiss$F_MISS > 0.6, "INDV"], 
            quote = FALSE, row.names = FALSE, 
            col.names = FALSE, "indvs_missing60.txt")

write.table(imiss[imiss$F_MISS > 0.5, "INDV"], 
            quote = FALSE, row.names = FALSE, 
            col.names = FALSE, "indvs_missing50.txt")

################################################################################
# SNP depth --------------------------------------------------------------------
################################################################################

SNPdepth <- read.table("out.ldepth.mean", header = T); colnames(SNPdepth)

(mean_snp_depth <- summary(SNPdepth$MEAN_DEPTH))
quantile(SNPdepth$MEAN_DEPTH, probs = c(0.1, 0.5, 0.9, 0.95, 0.99, 0.999, 0.9999))
nrow(SNPdepth[SNPdepth$MEAN_DEPTH > 100, ])

DEPTHLT <- 3
DEPTHUT <- 13

SNPdepth$filter <- case_when(SNPdepth$MEAN_DEPTH <  DEPTHLT | SNPdepth$MEAN_DEPTH >  DEPTHUT ~ "F",
                             SNPdepth$MEAN_DEPTH >= DEPTHLT | SNPdepth$MEAN_DEPTH <= DEPTHUT ~ "P")

(snp_depthhist <- ggplot() +
    geom_histogram(data = SNPdepth[SNPdepth$filter == "P",],
                   aes(MEAN_DEPTH), bins = 150, 
                   color = "black", fill = "gray80") +
    geom_histogram(data = SNPdepth[SNPdepth$filter == "F",],
                   aes(MEAN_DEPTH), bins = 150, 
                   color = "gray90", fill = "gray99") +
    scale_x_continuous(breaks = seq(0, 13, 2)) +
    theme_bw() + xlim(c(4,13)) +
    labs(x = "Mean SNP depth", y = "Frequency"))

ggsave("snp_depth.tiff", width = 12,
       height = 6, dpi = 300)

################################################################################
# SNP Missingness --------------------------------------------------------------
################################################################################

SNPmiss <- read.table("out.lmiss", header = T)

MISST <- 4/10
SNPmiss$filter <- case_when(SNPmiss$F_MISS >= MISST ~ "F",
                            SNPmiss$F_MISS <  MISST ~ "P")

(snp_miss_summ <- summary(SNPmiss$F_MISS))

(snp_misshist <- ggplot() +
    geom_histogram(data = SNPmiss[SNPmiss$filter == "P", ], 
                   aes(F_MISS), bins = 300, color = "black", 
                   fill = "gray80") +
    geom_histogram(data = SNPmiss[SNPmiss$filter == "F", ], 
                   aes(F_MISS), bins = 300, color = "gray90", 
                   fill = "gray99") +
    geom_vline(xintercept = MISST, color = "red") +
    labs(x = "SNP Missingness (%)", y = "Frequency") +
    theme_bw())


ggsave("snp_missingnessh.tiff", width = 12,
       height = 6, dpi = 300)


################################################################################
# SNPs Karyotype ---------------------------------------------------------------
################################################################################

# Read in genome organization stats (from NCBI: GCF_018296145.1).
genome <- read.table("../../info_files/sequence_report_Otsh_v2.0.tsv", 
                     sep = "\t", header = T) %>%
  # Retain autosomes only for present purposes.
  filter(Role != "unplaced-scaffold" & Chromosome.name != "MT") %>% 
  # Reorganize for reading into toGRanges (see ?toGRanges).
  select(c("Chromosome.name", "Seq.length")) %>% 
  `colnames<-`(., c("chr", "end")) %>% 
  # Relabel chromosomes to a consistent "Ots##" format.
  mutate(LG = paste0("Ots", gsub("^.{0,2}", "", chr)),
         start = 1, tpos = cumsum(end)) %>% 
  relocate(start, .after = chr)

# Plot base autosomes.
(base <- ggplot() +
    theme_minimal() +
    theme_light() +
    geom_bar(data = genome,
             aes(x = LG, y = end/1e6),
             stat = 'identity',
             fill = 'grey80',
             width = 1/5) +
    labs(x = NULL, y = "Base pair position (Mbp)") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          axis.ticks.x = element_blank(),
          panel.grid = element_blank()))

ggsave("chrom_snps.tiff", width = 15, height = 10, dpi = 300)


SNPs <- SNPdepth[,c("CHROM", "POS", "filter")] %>% 
  mutate(LG = case_when(
    CHROM == "NC_056429.1" ~ "Ots01", CHROM == "NC_056430.1" ~ "Ots02",
    CHROM == "NC_056431.1" ~ "Ots03", CHROM == "NC_056432.1" ~ "Ots04",
    CHROM == "NC_056433.1" ~ "Ots05", CHROM == "NC_056434.1" ~ "Ots06",
    CHROM == "NC_056435.1" ~ "Ots07", CHROM == "NC_056436.1" ~ "Ots08",
    CHROM == "NC_056437.1" ~ "Ots09", CHROM == "NC_056438.1" ~ "Ots10",
    CHROM == "NC_056439.1" ~ "Ots11", CHROM == "NC_056440.1" ~ "Ots12",
    CHROM == "NC_056441.1" ~ "Ots13", CHROM == "NC_056442.1" ~ "Ots14",
    CHROM == "NC_056443.1" ~ "Ots15", CHROM == "NC_056444.1" ~ "Ots16",
    CHROM == "NC_056445.1" ~ "Ots17", CHROM == "NC_056446.1" ~ "Ots18", 
    CHROM == "NC_056447.1" ~ "Ots19", CHROM == "NC_056448.1" ~ "Ots20",
    CHROM == "NC_056449.1" ~ "Ots21", CHROM == "NC_056450.1" ~ "Ots22",
    CHROM == "NC_056451.1" ~ "Ots23", CHROM == "NC_056452.1" ~ "Ots24",
    CHROM == "NC_056453.1" ~ "Ots25", CHROM == "NC_056454.1" ~ "Ots26",
    CHROM == "NC_056455.1" ~ "Ots27", CHROM == "NC_056456.1" ~ "Ots28",
    CHROM == "NC_056457.1" ~ "Ots29", CHROM == "NC_056458.1" ~ "Ots30",
    CHROM == "NC_056459.1" ~ "Ots31", CHROM == "NC_056460.1" ~ "Ots32",
    CHROM == "NC_056461.1" ~ "Ots33", CHROM == "NC_056462.1" ~ "Ots34"
  )) %>% 
  filter(!is.na(LG))

(snp_positions <- base +
  geom_point(data = SNPs, 
             aes(x = LG, y = POS/1e6), 
             color = "grey20",
             shape = 22, 
             size = 1/1000, 
             alpha = 2/5,
             position = position_jitter(1/5)))

ggsave("snp_positions.tiff", width = 15, height = 10, dpi = 300)


################################################################################
################################################################################