setwd("~/ots_landscape_genetics/pop_gen")

library(tidyverse); library(vcfR); library(ggplot2); library(adegenet)
library(viridis); library(poppr); library(hierfstat); library(R.utils)

nc <- 12 # Number of cores to use in parallel functions.

# Read in genetic data and convert to genlight object.
gunzip("../data/snps_maf2_m60.vcf.gz")
hdgbs <- read.vcfR("../data/snps_maf2_m60.vcf")
hdgl <- vcfR2genlight(hdgbs)

# Part I: Indv and pop info ----------------------------------------------------

sites <- read.delim(file = "../data/ch2023_sequenced.txt") %>% 
  arrange(Latitude) %>% 
  mutate(site = tools::toTitleCase(tolower(gsub("\\_.*", "", Population))),
         sitenum = as.numeric(rownames(.)))
sites[sites$site == "San", "site"] <- "San Juan"
sites[sites$site == "Salmon", "site"] <- "Salmon Fork"
sites[sites$site == "Big", "site"] <- "Big Salmon"

samples <- as.factor(c(colnames(hdgbs@gt[,c(2:ncol(hdgbs@gt))])))

# Read in sample and population info. Arrange accordingly. 
indvs <- read.csv("../data/landgen_chinook_indvs.csv", na.strings = "") %>% 
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

hdgl <- dartR::gl.sort(x = hdgl, order.by = c(sites$site))

# Part 2: Genetic diversity ----------------------------------------------------

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

#### Use pixy to est. nucleotide div. and diff. see https://pixy.readthedocs.io/en/latest/index.html ####
# Need WSL approval?

# Part 3a: Population structure via PCA ----------------------------------------

# hdgl_sub <- gl.drop.ind(hdgl)

# Takes a few hours.
hd_pca <- glPca(hdgl, nf = 3, parallel = T, n.cores = 14)

# Include glPCA object from earlier to decrease computation time.
# ncl <- find.clusters.genlight(hdgl, glPca = hd_sub_pca,
#                               max.n.clust = length(unique(hdgl@pop))/3)

write.csv(hd_pca$scores, "hdgbs_pcascores.csv", row.names = TRUE)

"%ni%" <- Negate("%in%")
drop_indvs <- (pc_scores[pc_scores$PC1 > 50, "fish_ID"])
hdgl_sub <- hdgl[indNames(hdgl) %ni% drop_indvs]
hd_sub_pca <- glPca(hdgl_sub, nf = 3, parallel = T, n.cores = 14)
write.csv(hd_sub_pca$scores, "hdgbs_sub_pcascores.csv", row.names = TRUE)

#pca_scores <- read.csv("../data/hdgbs_pcascores.csv")
pca_scores <- hd_sub_pca$scores

# Isolate eigenvalues (% variation explained for each PC axis).
(pc_var <- c(hd_pca$eig/sum(hd_pca$eig)*100)[1:5]); barplot(hd_pca$eig)

# Read in population coordinates for colours later on.
# Have to make some 3-letter codes into 4-letter codes. 
coords <- read.table("../data/ch2023_sequenced.txt", sep = "\t", header = T) %>% 
  mutate(short = str_to_title(substr(Population, start = 1, stop = 3))) 
coords[coords$Population == "TAKHANNE_RIVER", "short"] <- "Takha"
coords[coords$Population == "KILDALA_RIVER" , "short"] <- "Kild"


locs <- read.csv("../data/landgen_chinook_indvs.csv") %>% 
  mutate(site_full = gsub("[[:space:]]Brood", "", site_full))
locs[locs$site_full == "Upper Pitt", "site_full"] <- "Pitt"

# Isolate individual-based PCA vales and combine with geographical information.
pc_scores <- as.data.frame( hd_sub_pca$scores) %>% 
  rownames_to_column(var = "fish_ID") %>% 
  merge(., locs, by = "fish_ID") %>% 
  merge(., coords, by.x = "site_full", by.y = "Site") %>% 
  arrange(desc(Latitude)) %>% # For ggplot colouring. 
  mutate(site_full = factor(reorder(site_full, Latitude))) 

# Plot PCA scores and colour by population witcat ../01   ls -l
#h dark = south, light = north.
(snp_pca <- ggplot(data = pc_scores, 
       aes(x = PC1, y = PC2, group = site_full, fill = factor(Latitude))) +
  geom_point(shape = 21, size = 3/2) +
  # labs(x = paste0("PC1 (", round(pc_var[1], 1), "%)"),
  #      y = paste0("PC2 (", round(pc_var[2], 1), "%)")) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "right") +
  scale_fill_manual(values = alpha(c(viridis_pal(option = "D")(length(unique(pc_scores$site_full)))), 3/4),
                    labels = levels(pc_scores$site_full)) +
  guides(fill=guide_legend(ncol = 3, reverse = TRUE)))

ggsave("../plots/pca_fillscale_sub.tiff", dpi = 300, width = 8, height = 5)

ggplotly(snp_pca)


# Part 3b: Population structure via sNMF ---------------------------------------


# Heavily inspired by Connor French's 2020 population genetics guide:
# https://connor-french.github.io/intro-pop-structure-r/

library(LEA)

# gunzip("../data/hdgbs_snps_maf2_m60.vcf.gz")
LEA::vcf2geno(input.file  = "../data/hdgbs_snps_maf2_m60.vcf",
              output.file = "../data/hdgbs_snps_maf2_m60.geno")

maxK <- 10; rep <- 5

# Commented out! Only run this first time around, otherwise it overwrites and takes forever.
# ch_admix <- snmf(input.file = "../data/hdgbs_snps_maf2_m60.geno",
#                  entropy = TRUE, repetitions = rep, # Uses cross-entropy criterion for each repetition of K. 
#                  project = "new", K = 1:maxK,       # Use new project and test across 1-10 clusters.
#                  ploidy  = 2, CPU = nc,             # Use 12 cores and specify diploid genome.
#                  I = 3e2, seed = 240)               # Seed is my office number. Initialize with 3k SNPs.

# For future runs, just read in the project. 
ch_admix <- load.snmfProject("../data/hdgbs_snps_maf2_m60.snmfProject")

################ Adjust I above to something reasonable ########################
      ### Check defaults, use all SNPs?

plot(ch_admix, pch = 19)

# Generate an empty dataframe for populating.
snmf_vals <- as.matrix.data.frame(matrix(NA, ncol = maxK, nrow = rep)) %>% 
  `rownames<-`(., c(paste0("Run", seq(1, rep, 1)))) %>% 
  `colnames<-`(., c(paste0("K", seq(1, maxK, 1))))

# Populate above DF with CE values of each K value and each run. 
for (i in c(1:maxK)) {
  CE <- c(cross.entropy(ch_admix, K = i))
  snmf_vals[,i] <- CE
}

# Determine optimum run and K value.
opt_snmf <- which(snmf_vals == min(snmf_vals), arr.ind = T) %>% 
    as.data.frame() %>% `colnames<-`(., c("run", "K"))
paste("Optimum is K =", opt_snmf[2], "and run", opt_snmf[1])

# Organize q-value data into a organized dataframe.
qdf <- as.data.frame(LEA::Q(ch_admix, K = as.numeric(opt_snmf$K), 
                             run = as.numeric(opt_snmf$run))) %>% 
  mutate(fish = hdgl@ind.names) %>% relocate(fish, .before = "V1") %>% 
  `colnames<-`(., c("fish", paste0("K", seq(1:as.numeric(opt_snmf$K))))) %>% 
  pivot_longer(cols = starts_with("K"), names_to = "pop", values_to = "q") %>% 
  group_by(fish) %>% 
  mutate(likely_assignment = pop[which.max(q)],
         assignment_prob   = max(q)) %>% 
  arrange(likely_assignment, assignment_prob, q) %>% 
  ungroup() %>% 
  mutate(fish = forcats::fct_inorder(factor(fish)))



# https://teunbrand.github.io/ggh4x/reference/guide_axis_nested.html



(qcol <- ggplot(data = qdf) +
  geom_col(aes(x = fish, y = q, fill = pop), linewidth = 1) +
  theme(axis.text.x = element_blank()) +
  labs(x = NULL, y = "Ancestry coefficients") +
  theme_bw() +
  theme(panel.background = element_rect(color = "black", fill = NA),
        panel.spacing.x = unit(0, "lines"),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks.y = element_line(color = "black"),
        legend.margin = unit(c(0,0,0,0), "cm"),
        legend.position = "top",
        legend.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
  scale_y_continuous(expand = c(0,0)))

ggsave("../plots/q_cols.tiff", dpi = 300, width = 8, height = 5)

# Part 4: Genetic differentiation ----------------------------------------------


# Compute pairwise Fst values. Takes a while. Save locally.
gendist <- dartR::gl.fst.pop(hdgl, nclusters = nc)
write.csv(gendist, "../data/fst_matrix.csv", row.names = T)
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

ggsave("../plots/fst_heatmap.tiff", dpi = 300,
       width = 7, height = 7)

# FST summary.
summary(gdlf$Value)
hist(gdlf$Value, breaks = nrow(gendist), main = NULL, xlab = "FST")



# Part 5: NJTree ---------------------------------------------------------------

library(ape)
nj <- dartR::gl.tree.nj(x = hdgl)
plot(ape::unroot(nj), type = "unrooted", cex = 1, edge.width = 2, lab4ut = "axial")

