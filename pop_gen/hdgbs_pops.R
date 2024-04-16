setwd("~/ots_landscape_genetics/pop_gen")

library(tidyverse); library(vcfR); library(ggplot2); library(adegenet)
library(viridis); library(poppr); library(hierfstat)

# Read in genetic data and convert to genlight object.
hdgbs <- read.vcfR("../data/hdgbs_snps_maf2_m60.vcf.gz")
hdgl <- vcfR2genlight(hdgbs)


# Part I: Genetic diversity ----------------------------------------------------

# Isolate sample names.
samples <- as.factor(c(colnames(hdgbs@gt[,c(2:ncol(hdgbs@gt))])))

# Read in sample and population info. Arrange accordingly. 
indvs <- read.csv("../data/landgen_chinook_indvs.csv", na.strings = "") %>% 
  filter(fish_ID %in% c(colnames(hdgbs@gt[,c(2:ncol(hdgbs@gt))]))) %>% 
  arrange(factor(fish_ID, levels = samples))

# Fix one sample label for consistency (HarB == Har).
indvs[indvs$site_full == "Harrison Brood", "site_full"] <- "Harrison"
# ^ Are Harrison and Harrison Brood different? Doubt it but check w/ TH. 

# c(colnames(hdgbs@gt))[c(colnames(hdgbs@gt)) %ni% indvs$fish_ID]
# ^ Use this to check samples when pipeline is finished.

# Assign population factor.
hdgl@pop <- as.factor(indvs$site_full)

# Calculate genetic diversity summary stats.
snp_div <- dartR::gl.basic.stats(hdgl) # Takes a while.

# Isolate population level diversity stats.
pop_div <- data.frame(
  Ho  = round(apply(snp_div$Ho,  MARGIN = 2, FUN = mean, na.rm = T), 4),
  He  = round(apply(snp_div$Ho,  MARGIN = 2, FUN = mean, na.rm = T), 4),
  Fis = round(apply(snp_div$Fis, MARGIN = 2, FUN = mean, na.rm = T), 4)
) %>% rownames_to_column(var = "Sample")

rm(snp_div); gc() # Clear up memory
# per-SNP info takes up a lot of space that we don't need here. 

write.csv(pop_div, "../data/pop_snp_div.csv", row.names = FALSE)


# Part II: Population structure ------------------------------------------------


## NEEDS REDOING/REFINEMENT ONCE THE PIPELINE IS FINISHED ##



# hdgl_sub <- gl.drop.ind(hdgl)

# Takes a few hours.
# hd_pca <- glPca(hdgl, nf = 3, parallel = T, n.cores = 12)

# write.csv(hd_pca$scores, "hdgbs_pcascores.csv", row.names = TRUE)


"%ni%" <- Negate("%in%")
# drop_indvs <- rownames(pc_scores[pc_scores$PC1 > 50, ])
# hdgl_sub <- hdgl[indNames(hdgl) %ni% drop_indvs]
# hd_sub_pca <- glPca(hdgl_sub, nf = 3, parallel = T, n.cores = 12)
# write.csv(hd_sub_pca$scores, "hdgbs_sub_pcascores.csv", row.names = TRUE)

pca_scores <- read.csv("hdgbs_sub_pcascores.csv") %>% 
  column_to_rownames("X")

# Isolate eigenvalues (% variation explained for each PC axis).
(pc_var <- c(hd_sub_pca$eig/sum(hd_sub_pca$eig)*100)[1:5]); barplot(hd_pca$eig)

# Read in population coordinates for colours later on.
# Have to make some 3-letter codes into 4-letter codes. 
coords <- read.table("../map/ch2023_sequenced.txt", sep = "\t", header = T) %>% 
  mutate(short = str_to_title(substr(Population, start = 1, stop = 3))) %>% 
  filter(Population != "TAKHANNE_RIVER")
coords[coords$Population == "TAKHANNE_RIVER", "short"] <- "Takha"
coords[coords$Population == "KILDALA_RIVER" , "short"] <- "Kild"

# Isolate individual-based PCA vales and combine with geographical information.
pc_scores <- as.data.frame(pca_scores) %>% 
  mutate(pop = gsub("-.*","", rownames(.))) %>% 
  merge(., coords, by.x = "pop", by.y = "short") %>% 
  arrange(desc(Latitude)) %>% # For ggplot colouring. 
  mutate(pop = factor(reorder(pop, Latitude))) 

# Plot PCA scores and colour by population with dark = south, light = north.
ggplot(data = pc_scores, 
       aes(x = PC1, y = PC2,  fill = factor(Latitude), group = pop)) +
  geom_point(shape = 21, size = 3/2) +
  # labs(x = paste0("PC1 (", round(pc_var[1], 1), "%)"),
  #      y = paste0("PC2 (", round(pc_var[2], 1), "%)")) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "right") +
  scale_fill_manual(values = alpha(c(viridis_pal(option = "D")(length(unique(pc_scores$pop)))), 3/4),
                    labels = levels(pc_scores$pop)) +
  guides(fill=guide_legend(ncol = 3, reverse = TRUE))

ggsave("plots/pca_fillscale.tiff", dpi = 300, width = 8, height = 5)


# ------------------------------------------------------------------------------


