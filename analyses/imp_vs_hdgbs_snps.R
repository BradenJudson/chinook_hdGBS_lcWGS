setwd("~/ots_landscape_genetics/analyses/")

library(tidyverse); library(purrr); library(ggpmisc); library(ggpointdensity)

################################################################################
################# Part I: Imputed loci only ####################################
################################################################################

# Allele matching --------------------------------------------------------------

# Write a function for reading in and formatting frequency data per individual.
# File output from vcftools --freq
rf <- function(dir) paste0(dir, 
                    list.files(pattern = ".*.txt.frq", 
                    path = dir)) %>% set_names(.) %>% 
  map_df(read_table, .id = "FileName", 
         col_names = c("CHROM", 
                       "POS", "N_ALLELES", 
                       "N_CHR", "MajAll", 
                       "MinAll")) %>% 
  filter(MajAll != "{ALLELE:FREQ}") %>% 
  mutate(sample = gsub("_freq.txt.frq", "", gsub(dir, "", FileName)),
         AL1AF = as.numeric(sub(".*:", "", MajAll)),
         AL1Nuc  = as.factor(substr(MajAll, start = 0, stop = 1)),
         AL2AF = as.numeric(sub(".*:", "", MinAll)),
         AL2Nuc = as.factor(substr(MinAll, start = 0, stop = 1))) %>% 
  select(!c("MajAll", "MinAll", "N_CHR", "N_ALLELES", "FileName")) %>% 
  mutate(site = paste0(CHROM,"_",POS))

# Read in each dataset and format.
hdgbs <- rf(dir = "../data/ind_afs/hdgbs_n361_ref/") %>% 
  mutate(ref_frq_hd = AL1AF)
# hdgbs$ref_frq <- hdgbs$AL1AF
lcwgs <- rf(dir = "../data/ind_afs/impLC_n361_ref/") %>% 
  mutate(ref_frq_lcwg = AL1AF)
lcwgs$sample <- gsub("_freq_GLs.txt.frq", "", lcwgs$sample)
dat <- merge(hdgbs, lcwgs, by = c("CHROM", "POS", "sample", "site")) %>% 
  rename("ref_hd_nuc" = AL1Nuc.x, "alt_hd_nuc" = AL2Nuc.x,
         "ref_lc_nuc" = AL1Nuc.y, "alt_lc_nuc" = AL2Nuc.y,
         "ref_hd_frq" = AL1AF.x,  "alt_hd_frq" = AL2AF.x,
         "ref_lc_frq" = AL1AF.y,  "alt_lc_frq" = AL2AF.y)

# Individual-level dataframe that counts matched/mismatched alleles.
mismatches <- dat %>% group_by(sample) %>% 
  summarise(match     = sum(((ref_hd_nuc == ref_lc_nuc & alt_hd_nuc == alt_lc_nuc) & (ref_hd_frq == ref_lc_frq & alt_hd_frq == alt_lc_frq)) | ((ref_hd_nuc == ref_lc_nuc & alt_hd_nuc != alt_lc_nuc) & (ref_hd_frq == 1 & ref_lc_frq == 1))),
            n.match1  = sum(((ref_hd_nuc == ref_lc_nuc & alt_hd_nuc == alt_lc_nuc) & (abs(ref_hd_frq-ref_lc_frq) == 0.5 & abs(alt_hd_frq-alt_lc_frq) == 0.5)) | ((ref_hd_nuc == ref_lc_nuc & alt_hd_nuc != alt_lc_nuc) & (ref_hd_frq == 0.5 & ref_lc_frq == 0.5))),
            n.match2  = sum(((ref_hd_nuc == ref_lc_nuc & alt_hd_nuc == alt_lc_nuc) & (abs(ref_hd_frq-ref_lc_frq) == 1 & abs(alt_hd_frq-alt_lc_frq) == 1)) | ((ref_hd_nuc == ref_lc_nuc & alt_hd_nuc != alt_lc_nuc) & (ref_hd_frq == 0 & ref_lc_frq == 0)) | (ref_hd_nuc != ref_lc_nuc & alt_hd_nuc != alt_lc_nuc)),
            acc_match = round((match/n()*100), 2),
            acc_m1    = round((n.match1/n()*100), 2),
            acc_m2    = round((n.match2/n()*100), 2),
            total_check = round((match/n()*100) + (n.match1/n()*100) + (n.match2/n()*100), 1))


# Number of SNPs that match for both alleles, the reference allele but the alternate allele, etc.
length(unique(dat[(dat$ref_hd_nuc == dat$ref_lc_nuc) & (dat$alt_hd_nuc == dat$alt_lc_nuc), "site"]))
length(unique(dat[(dat$ref_hd_nuc == dat$ref_lc_nuc) & (dat$alt_hd_nuc != dat$alt_lc_nuc), "site"]))
length(unique(dat[(dat$ref_hd_nuc != dat$ref_lc_nuc) & (dat$alt_hd_nuc == dat$alt_lc_nuc), "site"]))
length(unique(dat[(dat$ref_hd_nuc != dat$ref_lc_nuc) & (dat$alt_hd_nuc != dat$alt_lc_nuc), "site"]))

# Isolate the locations of the SNPs where both alleles match among both datasets (>99.99%)
matched_alleles <- unique(dat[(dat$ref_hd_nuc == dat$ref_lc_nuc) & (dat$alt_hd_nuc == dat$alt_lc_nuc), "site"])
# write.csv(c.allele, "./data/allele_mismatch_imp.csv", row.names = F)

# Establish panel labels for facet wrap in the next step.
fwlabs <- c('acc_match' = "Matching genotypes",
  'acc_m1'    = "Genotypes with one mismatch",
  "acc_m2"    = "Genotypes with two mismatches")

# Plot individual  mismatches, etc.
(alLF <- mismatches[,c("sample", "acc_match", "acc_m1", "acc_m2")] %>% 
  pivot_longer(cols = c("acc_match", "acc_m1", "acc_m2")) %>% 
    mutate(name = factor(name, levels = c("acc_match", "acc_m1", "acc_m2"))) %>% 
  ggplot(aes(x = name, y = value/100)) + geom_violin() +
    geom_boxplot(outlier.alpha = 0, width = 0.4) +
    theme_bw() +
    facet_wrap(. ~ name, scales = "free", 
               labeller = as_labeller(fwlabs)) +
    labs(x = NULL, y = NULL) + 
    scale_y_continuous(labels = scales::percent, n.breaks = 3) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank()))
ggsave("../plots/imp_accuracy2.tiff", dpi = 300, width = 10, height = 6)

# Characterize matching/mismatching rate distributions
summary(mismatches$acc_match); hist(mismatches$acc_match, main = NULL, xlab = "Accuracy (%)")
summary(mismatches$acc_m1); hist(mismatches$acc_m1, main = NULL, xlab = "Accuracy (%)")
summary(mismatches$acc_m2); hist(mismatches$acc_m2, main = NULL, xlab = "Accuracy (%)")


# Allele frequencies -----------------------------------------------------------


# Write a function to rotate the .frq files from earlier.
# For each position, calculate the average major (specified) allele frequency
# and the number of alleles used in that calculation (i.e., n).
afcalc <- function(df, allele) df[,c("CHROM", "POS", "sample", allele)] %>% 
  pivot_wider(names_from = sample, values_from = allele) %>% 
  group_by(CHROM, POS) %>% 
  summarise(AF = mean(c_across(3:(ncol(.)-2)), na.rm = T),
            n  = sum(!is.na(c_across(3:(ncol(.)-2)))))
  
# Run of the hdGBS dataset focusing on the major allele.
hd_afs <- afcalc(df = dat[dat$site %in% matched_alleles,], 
                 allele = "ref_frq") 

# hd_pos <- paste0(hd_afs$CHROM, "_", hd_afs$POS)

# Same as above but rename the AF column to specify lc. 
lc_afs <- afcalc(df = dat[dat$site %in% matched_alleles,], 
                 allele = "ref_frq_lcwg") %>% 
  rename("AFlci" = AF)

# Merge allele frequency files together and reformat slightly.
# Only keep SNPs with 12 or more alleles present (>5 individuals).
# Based on the lowest number of individuals retained per population.
freqs <- merge(hd_afs, lc_afs, by = c("CHROM", "POS","n")) %>% 
  mutate(diff = abs(AF - AFlci)) %>% 
  filter(n > 5) %>% 
  mutate(site = paste0(CHROM, "_", POS)) %>% 
  filter(site %in% matched_alleles)

write.csv(freqs, "../data/allele_frequencies_n361_refalt.csv", row.names = F)
freqs <- read.csv("../data/allele_frequencies_n361_refalt.csv")

(frqs <- ggplot(data = freqs, aes(x = AF, y = AFlci)) +
  geom_pointdensity() +
  scale_colour_continuous(low = "gray80", high = "gray20") +
  theme_bw() + xlim(0,1) +
  scale_x_continuous(breaks = seq(0,1,0.25)) +
  scale_y_continuous(breaks = seq(0,1,0.25)) +
  theme(legend.position = "none") +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  stat_smooth(method = "lm", se = F, colour = "black") +
  stat_poly_eq(use_label(c("R2")), label.x = "right",
               label.y = "bottom") +
  labs(x = "hdGBS reference allele frequency",
       y = "Imputed lcWGS reference allele frequency")) 

ggsave("../plots/imputed_vs_hdgbs_maj_frqs.tiff", dpi = 300, height = 8, width = 8)

(frqlm <- lm(data = freqs, AFlci ~ AF)); summary(frqlm)
hist(frqlm$residuals, breaks = 200); summary(frqlm$residuals)
# q <- frqlm$residuals
# sum(q < 0)

(rhist <- ggplot(data = data.frame(res = frqlm$residuals),
                aes(x = res)) + theme_bw() +
  geom_density() +
  labs(x = "Residual", y = "Count"))

cowplot::plot_grid(plotlist = list(frqs, rhist), 
                   ncol = 1, labels = c("a)", "b)"),
                   rel_heights = c(3,1))

ggsave("../plots/test.tiff", dpi = 300,
       width = 10, height = 10)


# SNP Density ------------------------------------------------------------------

# Using vcftools --gzvcf *.vcf.gz --SNPdensity 1000

lcwgs_snpden <- read.delim("../data/chinook_lcwgs_snp_density.snpden", sep = "\t")
hdgbs_snpden <- read.delim("../data/hdgbs_snp_density.snpden", sep = "\t")

summary(lcwgs_snpden$VARIANTS.KB); hist(lcwgs_snpden$VARIANTS.KB)
summary(hdgbs_snpden$VARIANTS.KB); hist(hdgbs_snpden$VARIANTS.KB)
