# Set working directory.
setwd(dirname(rstudioapi::getSourceEditorContext()$path)); getwd()

library(tidyverse); library(vcfR); library(ggplot2); library(adegenet)
library(viridis); library(poppr); library(hierfstat); library(R.utils)

nc <- 12

# Read in genetic data and convert to genlight object.
gunzip("../09c_filter_branch2/snps_maf001.vcf.gz")
hdgbs <- read.vcfR("../09c_filter_branch2/snps_maf001.vcf")
hdgl <- vcfR2genlight(hdgbs)

# Part I: Indv and pop info ----------------------------------------------------

sites <- read.delim(file = "../01_info_files/ch2023_sequenced.txt") %>% 
  arrange(Latitude) %>% 
  mutate(site = tools::toTitleCase(tolower(gsub("\\_.*", "", Population))),
         sitenum = as.numeric(rownames(.)))
sites[sites$site == "San", "site"] <- "San Juan"
sites[sites$site == "Salmon", "site"] <- "Salmon Fork"
sites[sites$site == "Big", "site"] <- "Big Salmon"

samples <- as.factor(c(colnames(hdgbs@gt[,c(2:ncol(hdgbs@gt))])))

# Read in sample and population info. Arrange accordingly. 
indvs <- read.csv("../01_info_files/landgen_chinook_indvs.csv", na.strings = "") %>% 
  filter(fish_ID %in% c(colnames(hdgbs@gt[,c(2:ncol(hdgbs@gt))]))) %>% 
  arrange(factor(fish_ID, levels = samples))

# Rename some samples for consistent labeling downstream.
indvs[indvs$site_full == "Big Qualicum Brood", "site_full"] <- "Qualicum"
indvs[indvs$site_full == "Serpentine Brood", "site_full"] <- "Serpentine"
indvs[indvs$site_full == "Upper Pitt", "site_full"] <- "Pitt"
indvs[indvs$site_full == "Spius Brood", "site_full"] <- "Spius"
indvs[indvs$site_full == "Harrison Brood", "site_full"] <- "Harrison"
# ^ Are Harrison and Harrison Brood different? Doubt it but check w/ TH. 

# c(colnames(hdgbs@gt))[c(colnames(hdgbs@gt)) %ni% indvs$fish_ID]
# ^ Use this to check samples when pipeline is finished.

# Assign population factor.
hdgl@pop <- as.factor(indvs$site_full)

sites <- sites[sites$site %in% indvs$site_full, ]

hdgl <- dartR::gl.sort(x = hdgl, order.by = c(sites$site))


# Compute pairwise Fst values. Takes a while. Save locally.
gendist <- dartR::gl.fst.pop(hdgl_sub, nclusters = 16)
write.csv(gendist, "../01_info_files/fst_matrix.csv", row.names = T)
# gendist <- read.csv("../data/fst_matrix.csv", row.names = 1)