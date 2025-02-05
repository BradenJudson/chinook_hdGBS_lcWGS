setwd("~/ots_landscape_genetics/comp")

library(tidyverse); library(purrr); library(ggpmisc)


# To do ------------------------------------------------------------------------

# Compare genotype likelihoods with imputed vs. non-imputed.
# Do the "allele matching and AF comparisons with the entirety of the hdGBS and 
  # imputed datasets (i.e., don't go individual-by-individual). Boxplots!



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

freqs <- read.csv("./data/allele_frequencies.csv")[1:20000,]

library(ggpointdensity)

ggplot(data = freqs, aes(x = AFhd, y = AFlc)) +
  geom_pointdensity() +
  scale_colour_continuous(low = "gray80", high = "gray20") +
  theme_bw() +
  scale_x_continuous(breaks = seq(0,1,0.25)) +
  scale_y_continuous(breaks = seq(0,1,0.25)) +
  theme(legend.position = "none") +
  geom_smooth(method = "lm", se = FALSE, colour = "black")



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

# Deal with allele frequencies output from actual genotypes (via vcftools) first.
(frq_files <- paste0("./data/", list.files(pattern = ".*\\.frq", path = "./data"))[1:2])

files <- lapply(frq_files, read_table,
                col_names = c("CHROM", "POS", "N_ALLELES", 
                              "N_CHR", "MajAll", "MinAll")) %>% 
  `names<-`(.,c("hdGBS", "lcWGS_imputed")) %>% 
  lapply(., function(x) x %>% filter(MajAll != "{ALLELE:FREQ}") %>% 
           mutate(MajAF  = as.numeric(sub(".*:", "", MajAll)),
                  MajNuc = as.factor(substr(MajAll, start = 0, stop = 1)),
                  MinAF  = as.numeric(sub(".*:", "", MinAll)),
                  MinNuc = as.factor(substr(MinAll, start = 0, stop = 1)),
                  genPos = paste0(CHROM, "_", POS)) %>% 
           filter(genPos %in% hd_pos) %>% 
           select(!c("MajAll", "MinAll", "N_CHR", "N_ALLELES", "genPos")))


# Process allele frequencies derived from genotype likelihoods second.
# Formatting is different and requires some modification for consistency.
LCAFs <- read_table("../data/lcafs_maf005.mafs", skip = 1,
                           col_names = c("CHROM", "POS", "MinAF")) %>% 
  filter(!is.na(MinAF)) %>% mutate(MinAF = as.numeric(MinAF),
                                   MajAF = 1-MinAF,
                                   genPos = paste0(CHROM, "_", POS)) %>% 
  filter(genPos %in% hd_pos) %>% select(c("CHROM", "POS", "MajAF"))  

# Join allele frequency files together for all 3 datasets.
afs <- purrr::map_df(files, ~as.data.frame(.x), .id = "dataset") %>% 
  select(c(1:4)) %>% rbind(., LCAFs %>% mutate(dataset = "lcWGS")) %>% 
  pivot_wider(names_from = dataset, values_from = MajAF) 

write.csv(afs, "allele_frequencies.csv", row.names = F)

# Write a function plot allele frequencies for common SNPs.
af_comp <- \(x, y) {
  (plot <- ggplot(data = afs, aes(x = {{x}}, y = {{y}})) +
     theme_bw() + geom_point(alpha = 1/6)  +
     geom_smooth(method= "lm", colour = "deepskyblue3") +
     stat_poly_eq(use_label(c("R2")), label.x = "right", label.y = "bottom") +
    geom_abline(slope = 1, intercept = 0, colour = "deepskyblue3", linetype = 2) +
     labs(x = paste(gsub("_", " ", rlang::as_label(eval(parse(text=enquo(x)))[[2]])), "major allele frequency"),
          y = paste(gsub("_", " ", rlang::as_label(eval(parse(text=enquo(y)))[[2]])), "major allele frequency")))
  
    return(plot)
    
    }

# Join plots of interest together.
cowplot::plot_grid(plotlist = list(
  af_comp(x = hdGBS, y = lcWGS_imputed),
  af_comp(x = hdGBS, y = lcWGS),
  af_comp(x = lcWGS, y = lcWGS_imputed)), 
  ncol = 1, nrow = 3)

ggsave("../plots/afs_comparisons.png", dpi = 300, height = 16, width = 8)







# ################################################################################
# ################# Part II: All SNPs ############################################
# ################################################################################
# 
# # All SNPs ----------------------------------------------------------------
# 
# 
# rff <- function(x) x %>% filter(MajAll != "{ALLELE:FREQ}") %>% 
#   mutate(MajAF = as.numeric(sub(".*:", "", MajAll)),
#          MajNuc  = as.factor(substr(MajAll, start = 0, stop = 1)),
#          MinAF = as.numeric(sub(".*:", "", MinAll)),
#          MinNuc = as.factor(substr(MinAll, start = 0, stop = 1))) %>% 
#   select(!c("MajAll", "MinAll", "N_CHR", "N_ALLELES"))
# 
# 
# # LCo <- rff(x = read_table("../data/global_original_freqGLs.frq",
# #                   col_names = c("CHROM", "POS", "N_ALLELES", 
# #                                 "N_CHR", "MajAll", "MinAll")))
# 
# LCo <- read_table("../data/lcafs_maf005.mafs", skip = 1,
#                   col_names = c("CHROM", "POS", "MinAF")) %>% 
#   filter(!is.na(MinAF)) %>% mutate(MinAF = as.numeric(MinAF),
#                                    MajAF = 1-MinAF,
#                                    genPos = paste0(CHROM, "_", POS)) %>% 
#   filter(!rownames(.) %in% c(43830, 82753, 121301)) 
# 
# LCi <- rff(x = read_table("../data/global_imputed_freqGLs.frq",
#                           col_names = c("CHROM", "POS", "N_ALLELES", 
#                                         "N_CHR", "MajAll", "MinAll"))) %>% 
#   rename(iMajAF = "MajAF", iMajNuc = "MajNuc", 
#          iMinAF = "MinAF", iMinNuc = "MinNuc") %>% 
#   mutate(GP = paste0(CHROM, "_", POS))
# 
# 
# fAFs <- merge(LCo, LCi, by = c("CHROM", "POS")) %>% 
#   pivot_longer(cols = c(MajAF, iMajAF))  %>% 
#   mutate(dataset = case_when(
#     name == "MajAF"  ~ "lcWGS",
#     name == "iMajAF" ~ "Imputed lcWGS")) 
# 
# 
# (lc_lm <- ggplot(data = fAFs) +
#   geom_histogram(aes(x = value),  bins = 40,
#                 colour = "black", fill = "grey80") +
#   facet_wrap(~dataset, ncol = 1,
#              strip.position = "right") +
#   labs(x = "Major allele frequency",
#        y = "Number of SNPs") +
#   theme_bw() + 
#   theme(strip.background = element_blank()))
# 
# ggsave("../plots/majallele_frq_lowcov.tiff", dpi = 300, 
#        width = 10, height = 8)
# 
# (summ <- fAFs %>% 
#   group_by(dataset) %>% 
#   summarise(mean   = mean(value),
#             stdev  = sd(value),
#             median = median(value)))


