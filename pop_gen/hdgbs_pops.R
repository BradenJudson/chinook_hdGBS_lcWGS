setwd("~/ots_landscape_genetics/pop_gen")

library(tidyverse); library(vcfR); library(ggplot2); library(adegenet)
library(viridis); library(poppr); library(hierfstat)

nc <- 12 # Number of cores to use in parallel functions.

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

# Compute pairwise Fst values. Takes a while. Save locally.
gendist <- dartR::gl.fst.pop(hdgl, nclusters = nc)
write.csv(gendist, "../data/fst_matrix.csv", row.names = T)

# Convert to long formfor plotting purposes. 
gdlf <- gendist %>% 
  as.data.frame(row.names = c(rownames(gendist))) %>% 
  rownames_to_column("Var1") %>% 
  pivot_longer(-Var1, names_to = "Var2", values_to = "Value") %>% 
  mutate(Var1 = factor(Var1, levels = unique(rownames(gendist))), 
         Var2 = factor(Var2, levels = unique(rownames(gendist)))) %>% 
  filter(!is.na(Value))

# Plot pairwise FST heat map. 
ggplot(data = gdlf, aes(Var1, Var2)) +
  geom_tile(aes(fill = Value), colour = "grey80") +
  scale_fill_gradient(low = "skyblue1", high = "royalblue3") +
  theme_bw() + labs(x = NULL, y = NULL, 
                    fill = expression(paste(F[ST]))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = c(1/10, 17/20),
        legend.box.background = element_rect(colour = "black", fill = "white")) +
  guides(fill = guide_colorbar(frame.colour = "black", ticks.colour = "black"))

ggsave("../plots/fst_heatmap.tiff", dpi = 300,
       width = 7, height = 7)

# FST summary.
summary(gdlf$Value)
hist(gdlf$Value, breaks = nrow(gendist), main = NULL, xlab = "FST")



# NJTree -----------------------------------------------------------------------

library(ape)
nj <- dartR::gl.tree.nj(x = hdgl)
plot(ape::unroot(nj), type = "unrooted", cex = 1, edge.width = 2, lab4ut = "axial")


