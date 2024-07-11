setwd("~/ots_landscape_genetics/pop_gen")

library(tidyverse); library(vcfR); library(ggplot2); library(adegenet)
library(viridis); library(cowplot); library(hierfstat); library(R.utils)

nc <- 12 # Number of cores to use in parallel functions.

# Read in genetic data and convert to genlight object.
hdgbs <- read.vcfR("../data/snps_maf001_singletons.vcf")
hdgl <- vcfR2genlight(hdgbs)


# Individuals and sites --------------------------------------------------------

indvs <- read.csv("../data/landgen_chinook_indvs.csv", na.strings = "") %>% 
  mutate(site_full = gsub("[[:space:]]Brood", "", site_full))
indvs[indvs$site_full == "Upper Pitt", "site_full"] <- "Pitt"

sites <- read.delim(file = "../data/ch2023_sequenced.txt") %>% 
  mutate(site = tools::toTitleCase(tolower(gsub("\\_.*", "", Population))),
         sitenum = as.numeric(rownames(.))) %>% 
  arrange(Latitude) 

# PCA w/ plink -----------------------------------------------------------------

# First need to add PLINK to .Renviron, then run the PCA on the VCF. Much faster than alternatives.
system("plink.exe --vcf ../data/snps_maf001_singletons.vcf --aec --pca --out ../data/snps_maf001_singletons")

pc_scores <- read_table2("../data/snps_maf001_singletons.eigenvec", col_names = F) %>% 
  .[,2:ncol(.)] %>% # Remove unnecessary first column.
  `colnames<-`(., c("id", paste0("PC", 1:(ncol(.)-1)))) %>% 
  merge(., indvs, by.x = "id", by.y = "fish_ID") %>% 
  merge(., sites[,c("Site", "Latitude")], by.x = "site_full", 
           by.y = "Site", all.x = T) %>% 
  mutate(site_full = as.factor(site_full))

eigenval <- scan("../data/snps_maf001_singletons.eigenval")

pc_plot <- function(x, y) {
  
  xPC <- rlang::as_label(eval(parse(text=enquo(x)))[[2]])
  yPC <- rlang::as_label(eval(parse(text=enquo(y)))[[2]])
  
  pc_var <- eigenval/sum(eigenval)*100
  
  (pca <- ggplot(data = pc_scores,
                 aes(x = {{x}}, y = {{y}}, group = site_full, fill = factor(Latitude))) +
      geom_point(shape = 21, size = 3/2) + theme_bw() +
      theme(legend.title = element_blank(),
            legend.position = "none") +
      labs(x = paste0(xPC, " (", round(pc_var[as.numeric(str_sub(start = 3, end = 3, xPC))], 1), "%)"),
           y = paste0(yPC, " (", round(pc_var[as.numeric(str_sub(start = 3, end = 3, yPC))], 1), "%)")) +
      scale_fill_manual(values = c(viridis_pal(option = "D")(length(unique(pc_scores$site_full)))),
                        labels = levels(unique(pc_scores$site_full))))
  
  return(pca)
  
}

(snp_pca <- cowplot::plot_grid(plotlist = list(
  pc_plot(x = PC1, y = PC2),
  pc_plot(x = PC1, y = PC3),
  pc_plot(x = PC2, y = PC3) +
    theme(legend.position = "right")
), ncol = 3, rel_widths = c(1,1,1.9)))

ggsave("../plots/pca_fillscale_multi.tiff", dpi = 300, width = 16, height = 5)


star_plot <- function(x, y) {
  
  xPC <- rlang::as_label(eval(parse(text=enquo(x)))[[2]])
  yPC <- rlang::as_label(eval(parse(text=enquo(y)))[[2]])
  
  pc_var <- eigenval/sum(eigenval)*100
  
  pc_popav <- pc_scores %>% 
    group_by(site_full) %>% 
    summarise(V1A = mean(PC1),
              V2A = mean(PC2))
  
  pcavectors2 <- merge(pc_scores, pc_popav, by = "site_full") %>% 
    arrange(Latitude) 
  
  levels(pcavectors2$site_full) <- c(unique(pc_scores$site_full))
  
  
  (pca <- ggplot(data = pcavectors2, aes(colour = site_full)) +
      geom_segment(aes(x = V1A, y = V2A, xend = PC1, yend = PC2), linewidth = 3/4, alpha = 2/3) +
      scale_color_manual(values = c(viridis_pal(option = "D")(length(unique(pcavectors2$site_full)))),
                        labels = levels((pcavectors2$site_full))) +
      theme_bw() +
      # scale_color_gradient(guide = 'none') +
      ggrepel::geom_label_repel(data = pc_popav, max.overlaps = Inf, 
                                min.segment.length = 0, box.padding = 1/2,
                                aes(x = V1A, y = V2A, label = site_full), 
                                inherit.aes = FALSE) +
      labs(x = paste0(xPC, " (", round(pc_var[as.numeric(str_sub(start = 3, end = 3, xPC))], 1), "%)"),
           y = paste0(yPC, " (", round(pc_var[as.numeric(str_sub(start = 3, end = 3, yPC))], 1), "%)")) +
      theme(legend.position = "none", legend.title = element_blank(),
            legend.text = element_text(size = 11)))
  
  return(pca)
  
}

star_plot(x = PC1, y = PC2)

# Part I: Indv and pop info ----------------------------------------------------

# sites <- read.delim(file = "../data/ch2023_sequenced.txt") %>% 
#   arrange(Latitude) %>% 
#   mutate(site = tools::toTitleCase(tolower(gsub("\\_.*", "", Population))),
#          sitenum = as.numeric(rownames(.)))
# sites[sites$site == "San", "site"] <- "San Juan"
# sites[sites$site == "Salmon", "site"] <- "Salmon Fork"
# sites[sites$site == "Big", "site"] <- "Big Salmon"
# 
# samples <- as.factor(c(colnames(hdgbs@gt[,c(2:ncol(hdgbs@gt))])))
#  
# # Read in sample and population info. Arrange accordingly. 
# indvs <- read.csv("../data/landgen_chinook_indvs.csv", na.strings = "") %>% 
#   filter(fish_ID %in% c(colnames(hdgbs@gt[,c(2:ncol(hdgbs@gt))]))) %>% 
#   arrange(factor(fish_ID, levels = samples))
# 
# # Rename some samples for consistent labeling downstream.
# indvs[indvs$site_full == "Big Qualicum Brood", "site_full"] <- "Qualicum"
# indvs[indvs$site_full == "Serpentine Brood", "site_full"] <- "Serpentine"
# indvs[indvs$site_full == "Upper Pitt", "site_full"] <- "Pitt"
# indvs[indvs$site_full == "Spius Brood", "site_full"] <- "Spius"
# indvs[indvs$site_full == "Harrison Brood", "site_full"] <- "Harrison"
# # ^ Are Harrison and Harrison Brood different? Doubt it but check w/ TH. 
# 
# # c(colnames(hdgbs@gt))[c(colnames(hdgbs@gt)) %ni% indvs$fish_ID]
# # ^ Use this to check samples when pipeline is finished.
# 
# # Assign population factor.
# hdgl@pop <- as.factor(indvs$site_full)
# 
# sites <- sites[sites$site %in% indvs$site_full, ]
# 
# hdgl <- dartR::gl.sort(x = hdgl, order.by = c(sites$site))

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


# PCA Experimentation ----------------------------------------------------------

ggplot(data = pc_scores,
       aes(x = site_full, y = PC1)) +
  geom_point() +
  theme_bw() + labs(x = NULL) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


pc_popav <- pc_scores %>% 
  group_by(site_full) %>% 
  summarise(V1A = mean(PC1),
            V2A = mean(PC2))

pcavectors2 <- merge(pc_scores, pc_popav, by = "site_full")


ggplot(data = pcavectors2, aes(fill = Latitude, colour = Latitude)) +
  geom_segment(aes(x = V1A, y = V2A, xend = PC1, yend = PC2), linewidth = 3/4) +
  geom_point(aes(x = PC1, y = PC2),color = "black", shape = 21, size = 2) + theme_bw() +
  scale_color_gradient(guide = 'none') +
  ggrepel::geom_label_repel(data = pc_popav, max.overlaps = Inf, 
                            min.segment.length = 0, box.padding = 1/2,
                            aes(x = V1A, y = V2A, label = site_full), 
                            inherit.aes = FALSE) +
  labs(x = paste0("PC1 (", round(varPC1, 1), "%)"),
       y = paste0("PC2 (", round(varPC2, 1), "%)")) + 
  theme(legend.position = "top", legend.title = element_blank(),
        legend.text = element_text(size = 11))


pc_long <- pcavectors2 %>% 
  pivot_longer(cols = c("PC1", "PC2")) 

ggplot(data = pc_long, aes(x = site_full, y = value)) +
  geom_point(alpha = 2/3) + facet_wrap(~name, ncol = 1) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = NULL, y = "PC score")

ggplot(data = pc_popav, aes(x = V1A, y = V2A, label = site_full)) +
  geom_point() + theme_bw() +
    geom_label_repel(max.overlaps = Inf, 
                     min.segment.length = 0)



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

