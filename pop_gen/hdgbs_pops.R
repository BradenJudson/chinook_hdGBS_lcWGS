setwd("~/ots_landscape_genetics/pop_gen")

library(tidyverse); library(vcfR); library(ggplot2); library(adegenet)
library(viridis); library(cowplot); library(hierfstat); library(R.utils)

nc <- 12 # Number of cores to use in parallel functions.

# Read in genetic data and convert to genlight object.
hdgbs <- read.vcfR("../data/snps_maf001_singletons.vcf")
(hdgl <- vcfR2genlight(hdgbs))


# Individuals and sites --------------------------------------------------------

indvs <- read.csv("../data/landgen_chinook_indvs.csv", na.strings = "") %>% 
  mutate(site_full = gsub("[[:space:]]Brood", "", site_full))
indvs[indvs$site_full == "Upper Pitt", "site_full"] <- "Pitt"

sites <- read.delim(file = "../data/ch2023_sequenced.txt") %>% 
  mutate(site = tools::toTitleCase(tolower(gsub("\\_.*", "", Population))),
         sitenum = as.numeric(rownames(.))) %>% 
  arrange(Latitude) 


# Below seems to have a negligible impact.
# system("plink.exe --vcf ../data/snps_maf001_singletons.vcf --aec --indep-pairwise 500 5 0.2 --out ../data/snps_maf001LD")
# system("plink.exe --vcf ../data/snps_maf001_singletons.vcf --aec --extract ../data/snps_maf001LD.prune.in --pca --out ../data/snps_maf001_ldprune")


# Part 1: PCA w/ plink @ maf > 0.01 --------------------------------------------

# First need to add PLINK to .Renviron, then run the PCA on the VCF. Much faster than alternatives.
system("plink.exe --vcf ../data/snps_maf001_singletons.vcf --aec --pca 376 --out ../data/snps_maf001_singletons")

format_eigenvec <- \(eigenvec_file) {
  df <- read.table(eigenvec_file) %>% .[,2:ncol(.)] %>%         # Remove unnecessary first column.
    `colnames<-`(., c("id", paste0("PC", 1:(ncol(.)-1)))) %>%   # Rename columns.
    merge(., indvs, by.x = "id", by.y = "fish_ID") %>%          # Add populations.
    merge(., sites[,c("Site", "Latitude")], by.x = "site_full", # Add site latitudes.
          by.y = "Site", all.x = T) %>%                         # Keep all individuals.
    mutate(site_full = as.factor(site_full))                    # For plotting later.
  return(df)
}

# Read in corresponding PC values and attach population-level information.
pc_scores001 <- format_eigenvec("../data/snps_maf001_singletons.eigenvec")

# Eigenvalues (for %var explained by each axis).
eigenval001 <- scan("../data/snps_maf001_singletons.eigenval")

# Function for plotting PC scores.
pc_plot <- function(maf, x, y) {
  
  # Isolate variables as either strings or objects as needed.
  xPC <- rlang::as_label(eval(parse(text=enquo(x)))[[2]])
  yPC <- rlang::as_label(eval(parse(text=enquo(y)))[[2]]) 
  scores_df <- get(paste0("pc_scores", maf), envir = .GlobalEnv)
  eigenvals <- get(paste0("eigenval",  maf), envir = .GlobalEnv)
  
  # Percent of total variation explained by each axis (approximate).
  pc_var <- eigenvals/sum(eigenvals)*100
  
  # Construct the plot itself.
  (pca <- ggplot(data = scores_df,
                 aes(x = {{x}}, y = {{y}}, group = site_full, 
                     fill = factor(Latitude))) +
      geom_point(shape = 21, size = 3/2) + theme_bw() +
      theme(legend.title = element_blank(),
            legend.position = "bottom") + guides(fill = guide_legend(nrow = 5, byrow = TRUE)) +
      labs(x = paste0(xPC, " (", round(pc_var[as.numeric(str_sub(start = 3, end = 3, xPC))], 2), "%)"),
           y = paste0(yPC, " (", round(pc_var[as.numeric(str_sub(start = 3, end = 3, yPC))], 2), "%)")) +
      scale_fill_manual(values = c(viridis_pal(option = "D")(length(unique(scores_df$site_full)))),
                        labels = levels(unique(scores_df$site_full))))
  return(pca)
  }

# Return all combinations of PCs 1-3.
(pca005 <- ggpubr::ggarrange( pc_plot(maf = "001", x = PC1, y = PC2),
                              pc_plot(maf = "001", x = PC1, y = PC3),
                              pc_plot(maf = "001", x = PC2, y = PC3),
                              align = "h", ncol = 3, legend = "top",
                              common.legend = TRUE, labels = c("a)")))


ggsave("../plots/hdgbs_pca_maf001.tiff", dpi = 300, width = 12, height = 6)


# Part 2: PCA w/ plink @ maf > 0.05 --------------------------------------------

system("plink.exe --vcf ../data/snps_maf001_singletons.vcf --maf 0.05 --aec --pca 376 --out ../data/snps_maf005_singletons")

pc_scores005 <- format_eigenvec("../data/snps_maf005_singletons.eigenvec")

eigenval005 <- scan("../data/snps_maf005_singletons.eigenval")

(pca001 <- ggpubr::ggarrange( pc_plot(maf = "005", x = PC1, y = PC2),
                              pc_plot(maf = "005", x = PC1, y = PC3),
                              pc_plot(maf = "005", x = PC2, y = PC3),
                              align = "h", ncol = 3, legend = "none",
                              labels = c("b)")))

ggsave("../plots/hdgbs_pca_maf005.tiff", dpi = 300, width = 16, height = 5)

cowplot::plot_grid(pca005, pca001, ncol = 1, rel_heights = c(1.4, 1))
ggsave("../plots/hdgbs_pca_full.tiff", dpi = 300, width = 12, height = 10)

# Part 2: Genetic diversity ----------------------------------------------------
# 
# # Calculate genetic diversity summary stats.
# snp_div <- dartR::gl.basic.stats(hdgl) # Takes a while.
# 
# # Isolate population level diversity stats.
# pop_div <- data.frame(
#   Ho  = round(apply(snp_div$Ho,  MARGIN = 2, FUN = mean, na.rm = T), 4),
#   He  = round(apply(snp_div$Ho,  MARGIN = 2, FUN = mean, na.rm = T), 4),
#   Fis = round(apply(snp_div$Fis, MARGIN = 2, FUN = mean, na.rm = T), 4)
# ) %>% rownames_to_column(var = "Sample")
# 
# rm(snp_div); gc() # Clear up memory
# # per-SNP info takes up a lot of space that we don't need here. 
# 
# write.csv(pop_div, "../data/pop_snp_div.csv", row.names = FALSE)
# 
# #### Use pixy to est. nucleotide div. and diff. see https://pixy.readthedocs.io/en/latest/index.html ####
# # Need WSL approval?