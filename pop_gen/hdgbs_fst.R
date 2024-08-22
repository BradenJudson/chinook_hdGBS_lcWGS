#!/usr/bin/env Rscript

for (p in c("tidyverse", "vcfR", "dartR")) {
  if (!suppressMessages(require(p, character.only = T))) {
    message(paste("Installing:", p))
    install.packages(p, repos = "https://mirror.its.dal.ca/cran", dependencies = T)}
}

hdgbs <- read.vcfR("../09b_filter_branch1/snps_maf001.vcf")
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


# c(colnames(hdgbs@gt))[c(colnames(hdgbs@gt)) %ni% indvs$fish_ID]
# ^ Use this to check samples when pipeline is finished.

# Assign population factor.
hdgl@pop <- as.factor(indvs$site_full)

sites <- sites[sites$site %in% indvs$site_full, ]


# Part 2: Genetic differentiation ----------------------------------------------

# Compute pairwise Fst values. Takes a while. Save locally.
gendist <- dartR::gl.fst.pop(hdgl_sub, nclusters = 16)
write.csv(gendist, "../12_fst/fst_matrix.csv", row.names = T)
# gendist <- read.csv("../data/fst_matrix.csv", row.names = 1)


# Convert to long formfor plotting purposes. 
gdlf <- gendist %>% 
  as.data.frame(row.names = c(rownames(.))) %>% 
  rownames_to_column("Var1") %>% 
  pivot_longer(-Var1, names_to = "Var2", values_to = "Value") %>% 
  mutate(Var1 = factor(Var1, levels = rownames(gendist)), 
         Var2 = factor(Var2, levels = rownames(gendist))) %>% 
  filter(!is.na(Value) & !is.na(Var2))

# Plot pairwise FST heat map. 
ggplot(data = gdlf, aes(Var1, Var2)) +
  geom_tile(aes(fill = Value), colour = "grey80") +
  scale_fill_gradient(low = "#12C7EB5E", high = "#0F37D6D8") +
  theme_bw() + labs(x = NULL, y = NULL, 
                    fill = expression(paste(F[ST]))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1/5),
        axis.text   = element_text(size = 7),
        legend.position = c(1/10, 17/20),
        legend.box.background = element_rect(colour = "black", fill = "white"),
        panel.grid = element_line(colour = "gray97")) +
  guides(fill = guide_colorbar(frame.colour = "black", ticks.colour = "black"))

ggsave("../12_fst/fst_heatmap.tiff", dpi = 300,
       width = 7, height = 7)

# FST summary.
summary(gdlf$Value)
hist(gdlf$Value, breaks = nrow(gendist), main = NULL, xlab = "FST")

