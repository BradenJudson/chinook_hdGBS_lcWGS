setwd("~/ots_landscape_genetics/hdGBS")

library(tidyverse); library(ggplot2); library(vcfR)

# Load HDplot function from https://github.com/gjmckinney/HDplot
# Also see McKinney et al., 2017 (DOI: 10.1111/1755-0998.12613).
source("./HDplot/HDplot.R")

input <- read.vcfR("../data/snps_maf001.vcf")

# Tresholds in Chinook described in McKinney et al., 2017.
Hmax <- 0.6
Dmax <- 7.0

# Add a qualifier column for paralogs or non-paralogs.
HDplotResults <- HDplot(input) %>% 
  mutate(SNP = as.factor(case_when(H >  0.60 | abs(D) >  7 ~ "Duplicate",
                                   H <= 0.60 | abs(D) <= 7 ~ "Singleton")))

summary(HDplotResults$SNP); head(HDplotResults)

(HD_HD <- ggplot(data = HDplotResults) +
    geom_vline(xintercept = Hmax, linewidth = 1,
               color = "gray70") +
    geom_hline(yintercept = c(Dmax, Dmax*-1), 
               linewidth = 1, color = "gray70") +
  geom_point(aes(x = H, y = D, fill = SNP, 
                 color = SNP), 
                 alpha = 2/5) +
    scale_color_manual(values = c("red3", "blue3")) +
    labs(x = "H", y = "D (read-ratio deviation)") +
    theme_bw() + theme(legend.position = "none",
    panel.grid = element_blank()))

ggsave("./stats/hdplot_fig.tiff", dpi = 300,
       width = 10, height = 8)

# Write putative paralog positions to directory. 
# Have to format as per PLINK requirements. Column 4 is arbitrary.
write.table(HDplotResults[HDplotResults$SNP == "Duplicate", c("CHROM", "POS")] %>% 
              mutate(POS2 = POS, AR = rownames(.)),
            "./stats/duplicateSNP_IDs.txt", sep = "\t", quote = F, 
            row.names = F, col.names = F)

# Filter using above.
system("plink.exe --vcf ../data/snps_maf001.vcf --aec --exclude range ./stats/duplicateSNP_IDs.txt --recode vcf --out ../data/snps_maf001_singletons")
