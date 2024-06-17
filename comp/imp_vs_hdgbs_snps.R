setwd("~/ots_landscape_genetics/comp")

library(tidyverse); library(purrr); library(ggpmisc)


# To do ------------------------------------------------------------------------

# Compare genotype likelihoods with imputed vs. non-imputed.
# Do the "allele matching and AF comparisons with the entirety of the hdGBS and 
  # imputed datasets (i.e., don't go individual-by-individual). Boxplots!



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
         MajAF = as.numeric(sub(".*:", "", MajAll)),
         MajNuc  = as.factor(substr(MajAll, start = 0, stop = 1)),
         MinAF = as.numeric(sub(".*:", "", MinAll)),
         MinNuc = as.factor(substr(MinAll, start = 0, stop = 1))) %>% 
  select(!c("MajAll", "MinAll", "N_CHR", "N_ALLELES", "FileName"))

# Read in hdGBS data.
hdgbs <- rf(dir = "./hdgbs_sub_frq/") 
  
# Read in lcWGS data. Renaming helps next step.
lcwgs <- rf(dir = "./imputed_snp_frq/") %>% 
  rename(MajAFlc = "MajAF", MajNuclc = "MajNuc", 
         MinAFlc = "MinAF", MinNuclc = "MinNuc")

write.table(lcwgs, quote = F, row.names = FALSE, "lcwgs_imp.txt")

# sum(hdgbs$RefNuc == lcwgs$RefNuclc)/nrow(lcwgs)*100

# Combine lcWGS and hdGBS datasets. 
dat <- merge(hdgbs, lcwgs, by = c("CHROM", "POS", "sample"))


# Rearrange columns such that major and minor alleles are consistent.
dat[dat$MajNuc==dat$MinNuclc & dat$MinNuc==dat$MajNuclc,  ] <- dat[dat$MajNuc==dat$MinNuclc & dat$MinNuc==dat$MajNuclc,   c(1:3,6:7,4:5,8:11) ]
dat[(dat$MajNuc==dat$MinNuclc & dat$MinNuc!=dat$MajNuclc),] <- dat[(dat$MajNuc==dat$MinNuclc & dat$MinNuc!=dat$MajNuclc), c(1:3,4:7,10:11,8:9)]
dat[(dat$MajNuc!=dat$MinNuclc & dat$MinNuc==dat$MajNuclc),] <- dat[(dat$MajNuc!=dat$MinNuclc & dat$MinNuc==dat$MajNuclc), c(1:3,6:7,4:5,8:11) ]
dat[(dat$MajNuc!=dat$MajNuclc & dat$MinNuc==dat$MinNuclc),] <- dat[(dat$MajNuc!=dat$MajNuclc & dat$MinNuc==dat$MinNuclc), c(1:3,6:7,4:5,10:11,8:9)]


# confirm same base pairs
(sb2 <- sum(dat$MajNuc == dat$MajNuclc & dat$MinNuc == dat$MinNuclc))
(sb1 <- sum(dat$MajNuc == dat$MajNuclc & dat$MinNuc != dat$MinNuclc))
(sb0 <- sum(dat$MajNuc != dat$MajNuclc & dat$MinNuc != dat$MinNuclc))
nrow(dat); (sb0 + sb1 + sb2); nrow(dat) - (sb0 + sb1 + sb2) == 0

# Allele matches, mismatches and accuracy by individual. 
c.allele <- dat %>% group_by(sample) %>% 
  summarise(match     = sum(((MajNuc == MajNuclc & MinNuc == MinNuclc) & (MajAF == MajAFlc & MinAF == MinAFlc)) | ((MajNuc == MajNuclc & MinNuc != MinNuclc) & (MajAF == 1 & MajAFlc == 1))),
            n.match1  = sum(((MajNuc == MajNuclc & MinNuc == MinNuclc) & (abs(MajAF-MajAFlc) == 0.5 & abs(MinAF-MinAFlc) == 0.5)) | ((MajNuc == MajNuclc & MinNuc != MinNuclc) & (MajAF == 0.5 & MajAFlc == 0.5))),
            n.match2  = sum(((MajNuc == MajNuclc & MinNuc == MinNuclc) & (abs(MajAF-MajAFlc) == 1 & abs(MinAF-MinAFlc) == 1)) | ((MajNuc == MajNuclc & MinNuc != MinNuclc) & (MajAF == 0 & MajAFlc == 0)) | (MajNuc != MajNuclc & MinNuc != MinNuclc)),
            acc_match = round((match/n()*100), 2),
            acc_m1    = round((n.match1/n()*100), 2),
            acc_m2    = round((n.match2/n()*100), 2),
            total_check = round((match/n()*100) + (n.match1/n()*100) + (n.match2/n()*100), 1))

write.csv(c.allele, "./data/allele_mismatch_imp.csv", row.names = F)

# Establish panel labels for facet wrap in the next step.
fwlabs <- c('acc_match' = "Percentage of matching genotypes",
  'acc_m1'    = "Percentage of genotypes with one mismatch",
  "acc_m2"    = "Percentage of genotypes with two mismatches")

(alLF <- c.allele[,c("sample", "acc_match", "acc_m1", "acc_m2")] %>% 
  pivot_longer(cols = c("acc_match", "acc_m1", "acc_m2")) %>% 
    mutate(name = factor(name, levels = c("acc_match", "acc_m1", "acc_m2"))) %>% 
  ggplot(aes(x = name, y = value/100)) + geom_boxplot(outlier.alpha = 0) +
    geom_jitter(width = 1/8, shape = 21,
                colour = "black", fill = "grey25", 
                alpha = 1/4) + theme_bw() +
    facet_wrap(. ~ name, scales = "free", labeller = as_labeller(fwlabs)) +
    labs(x = NULL, y = NULL) + 
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank()))

ggsave("../plots/imp_accuracy.tiff", dpi = 300, width = 10, height = 6)

summary(c.allele$acc_match); hist(c.allele$acc_match, main = NULL, xlab = "Accuracy (%)")
summary(c.allele$acc_m1); hist(c.allele$acc_m1, main = NULL, xlab = "Accuracy (%)")
summary(c.allele$acc_m2); hist(c.allele$acc_m2, main = NULL, xlab = "Accuracy (%)")

# Non-individual / global summary.
tot.acc <- dat %>% 
  summarise(match     = sum(((MajNuc == MajNuclc & MinNuc == MinNuclc) & (MajAF == MajAFlc & MinAF == MinAFlc)) | ((MajNuc == MajNuclc & MinNuc != MinNuclc) & (MajAF == 1 & MajAFlc == 1))),
            n.match1  = sum(((MajNuc == MajNuclc & MinNuc == MinNuclc) & (abs(MajAF-MajAFlc) == 0.5 & abs(MinAF-MinAFlc) == 0.5)) | ((MajNuc == MajNuclc & MinNuc != MinNuclc) & (MajAF == 0.5 & MajAFlc == 0.5))),
            n.match2  = sum(((MajNuc == MajNuclc & MinNuc == MinNuclc) & (abs(MajAF-MajAFlc) == 1 & abs(MinAF-MinAFlc) == 1)) | ((MajNuc == MajNuclc & MinNuc != MinNuclc) & (MajAF == 0 & MajAFlc == 0)) | (MajNuc != MajNuclc & MinNuc != MinNuclc)),
            acc_match = round((match/n()*100), 2),
            acc_m1    = round((n.match1/n()*100), 2),
            acc_m2    = round((n.match2/n()*100), 2),
            total_check = round((match/n()*100) + (n.match1/n()*100) + (n.match2/n()*100), 1))


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
hd_afs <- afcalc(df = dat, allele = "MajAF") 

hd_pos <- paste0(hd_afs$CHROM, "_", hd_afs$POS)

# Same as above but rename the AF column to specify lc. 
lc_afs <- afcalc(df = dat, allele = "MajAFlc") %>% 
  rename("AFlc" = AF)

# Merge allele frequency files together and reformat slightly.
freqs <- merge(hd_afs, lc_afs, by = c("CHROM", "POS")) %>% 
  rename("AFhd" = AF, "n_samples" = n.x) %>% select(c(1:5)) %>% 
  mutate(diff = abs(AFhd - AFlc)) %>% 
  filter(n_samples > 10) # To minimize bias of one-offs and rare alleles.

write.csv(freqs, "./data/allele_frequencies.csv", row.names = F)

(lcimp_hdgbs.lm <- ggplot(data = freqs, aes(x = AFhd, y = AFlc)) + 
  geom_point(alpha = 1/4) + theme_bw() + 
  labs(x = "hdGBS major allele frequency",
       y = "lcWGS major allele frequency") +
  geom_smooth(method = "lm", linetype = 2) +
  stat_poly_eq(use_label(c("R2")),label.x = "right",
                         label.y = "bottom"))

ggsave("../plots/imputed_vs_hdgbs_maj_frqs.tiff", dpi = 300, height = 8, width = 12)

(frqlm <- lm(data = freqs, AFlc ~ AFhd)); summary(frqlm)

hist(freqs$diff, main = NULL, xlab = "Difference in allele frequency")


# AF Comparisons ---------------------------------------------------------------

(frq_files <- paste0("./data/", list.files(pattern = ".*\\.frq", path = "./data")))

files <- lapply(frq_files, read_table,
                col_names = c("CHROM", "POS", "N_ALLELES", 
                              "N_CHR", "MajAll", "MinAll")) %>% 
  `names<-`(.,c("hdGBS", "lcWGS_imputed", "lcWGS")) %>% 
  lapply(., function(x) x %>% filter(MajAll != "{ALLELE:FREQ}") %>% 
           mutate(MajAF  = as.numeric(sub(".*:", "", MajAll)),
                  MajNuc = as.factor(substr(MajAll, start = 0, stop = 1)),
                  MinAF  = as.numeric(sub(".*:", "", MinAll)),
                  MinNuc = as.factor(substr(MinAll, start = 0, stop = 1)),
                  genPos = paste0(CHROM, "_", POS)) %>% 
           filter(genPos %in% hd_pos) %>% 
           select(!c("MajAll", "MinAll", "N_CHR", "N_ALLELES", "genPos")))

afs <- purrr::map_df(files, ~as.data.frame(.x), .id = "dataset") %>% 
  select(c(1:4)) %>% 
  pivot_wider(names_from = dataset, values_from = MajAF) %>% 
  pivot_longer(cols = c("lcWGS", "hdGBS")) %>% 
  mutate(diff = abs(lcWGS_imputed - value))

write.csv(afs, "allele_frequencies.csv", row.names = F)

(lm1 <- lm(data = afs[afs$name == "hdGBS",], lcWGS_imputed ~ value)); summary(lm1)

(afp <- ggplot(data = afs, aes(x = lcWGS_imputed, y = value, group = name)) +
  theme_bw() +stat_poly_eq(use_label(c("R2")),label.x = "right", label.y = "bottom") +
  geom_point(alpha = 1/6) +
  facet_grid(~ name) +
  geom_abline(slope = 1, intercept = 0, colour = "deepskyblue3", linetype = 2) +
  geom_smooth(method= "lm", colour = "deepskyblue3") +
  labs(x = "Imputed lcWGS major allele frequency",
       y = "Major allele frequency"))

ggsave("../plots/imputation_full.tiff", dpi = 300, height = 8, width = 15)

(lc_hd <- purrr::map_df(files[c("hdGBS", "lcWGS")], 
                       ~as.data.frame(.x), .id = "dataset") %>%
  select(c(1:4)) %>% 
  pivot_wider(names_from = dataset, values_from = MajAF) %>% 
  ggplot(data = ., aes(x = lcWGS, y = hdGBS)) + 
    geom_point(alpha = 1/6) +
    geom_abline(slope = 1, intercept = 0, colour = "deepskyblue3", linetype = 2) +
  theme_bw() + geom_smooth(method = "lm", se = FALSE,
                           colour = "deepskyblue3") +
  labs(x = "lcWGS major allele frequency",
       y = "hdGBS major allele frequency") +
  stat_poly_eq(use_label(c("R2")),label.x = "right",
               label.y = "bottom"))

ggsave("../plots/hdGBS.tiff", dpi = 300, height = 8, width = 12)

cowplot::plot_grid(lc_hd, afp, ncol = 1) +
  annotate(geom = 'text', x = 0.9, y = 0.065, 
           label = bquote("R^2 == 0.99"), parse = TRUE)

ggsave("../plots/afs_comparisons.tiff", dpi = 300, height = 12, width = 12)
