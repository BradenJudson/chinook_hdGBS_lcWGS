setwd("~/ots_landscape_genetics/lcWGS/imputation")

library(tidyverse)

hdgbs <- read.table("eqHDGBS.frq", header = T, row.names = NULL) %>% 
  `colnames<-`(., c("CHROM", "POS", "NALLELES", "N_CHR", "MajAll", "MinAll")) %>% 
  mutate(MAF = as.numeric(sub(".*:", "", MinAll)))

impv <- read.table("eqIMPVCF.frq", header = T, row.names = NULL) %>% 
  `colnames<-`(., c("CHROM", "POS", "NALLELES", "N_CHR", "MajAll", "MinAll")) %>% 
  mutate(MAF = as.numeric(sub(".*:", "", MinAll)))

comb <- merge(hdgbs[,c(1:2, 7)] %>% `colnames<-`(., c("CHROM", "POS", "HD_MAF")),
              impv[,c(1:2,  7)] %>% `colnames<-`(., c("CHROM", "POS", "IM_MAF")), 
              by = c("CHROM", "POS")) %>% 
  mutate(diff = abs(HD_MAF - IM_MAF))
(lm1 <- lm(data = comb, HD_MAF ~ IM_MAF)); summary(lm1)

summary(comb$diff); hist(comb$diff, xlab = "Difference in MAF")

ggplot(data = comb, aes(x = diff)) + 
  geom_histogram(colour = "black", fill = "grey90", bins = 100) +
  theme_bw() +
  labs(x = "Δ Minor allele frequency [abs(lcWGS - hdGBS)]",
       y = "Frequency")

ggplot(data = comb, aes(x = HD_MAF, y = IM_MAF)) + 
  geom_point(alpha = 1/4) + theme_bw() +
  geom_smooth(method = "lm", linetype = 2, colour = "red3") +
  labs(x = "Minor allele frequency (hdGBS)",
       y = "Minor allele frequency (imputed lcWGS)") +
  scale_x_continuous(breaks = seq(0, 0.5, 0.1)) +
  scale_y_continuous(expand = expand_scale(mult = c(0, 1/10))) +
  stat_poly_eq(use_label(c("R2", "p")),
               small.p = TRUE)

ggsave("../../plots/imputation_test.tiff", dpi = 300, width = 12, height = 10)

diff_file <- read.table("diffvcf.diff.sites_in_files", header = T) %>% 
  mutate(match = as.factor(case_when(ALT1 == ALT2 ~ "Y",
                                     ALT1 != ALT2 ~ "N")))

(sdf <- summary(diff_file$match))
sdf[1]/(sdf[1] + sdf[2])*100

