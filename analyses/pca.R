setwd("~/ots_landscape_genetics/analyses")

library(tidyverse); library(viridis); library(cowplot)

indvs <- read.csv("../data/landgen_chinook_indvs.csv", na.strings = "") %>% 
  mutate(site_full = gsub("[[:space:]]Brood", "", site_full))
indvs[indvs$site_full == "Upper Pitt", "site_full"] <- "Pitt"

# Retain populations shared by all datasets.
shared_pops <- unique(read.csv("../data/shared_samples_n361.csv")[,"site_full"]) 
sites <- read.delim(file = "../data/ch2023_sequenced.txt") %>% 
  mutate(Site = gsub("[[:space:]]", "", Site),
         sitenum = as.numeric(rownames(.))) %>% 
  arrange(Latitude) %>% 
  filter(Site %in% shared_pops)

# --pca 1000 just ensures that all possible PCs are calculated (so %var explained is accurate).
system("plink.exe --vcf ../data/vcfs_n361/hdgbs_maf5_m70_pruned.vcf.gz  --aec --double-id --pca 1000 --out ../data/pca/hdgbs_pruned")
system("plink.exe --vcf ../data/vcfs_n361/lcwgs_ldpruned_maf005_n361.vcf.gz --aec --double-id --pca 1000 --out ../data/pca/lcwgs_imputed_pruned")

# Format PCA output data from above.
format_eigenvec <- \(eigenvec_file) {
  df <- read.table(eigenvec_file) %>% .[,2:ncol(.)] %>%            # Remove unnecessary first column.
    `colnames<-`(., c("id", paste0("PC", 1:(ncol(.)-1)))) %>%      # Rename columns.
    mutate(id = gsub("309\\-", "309\\.", gsub('^.*\\_[0-9.]+', "", # Make label formatting consistent.
                gsub(pattern = "\\.dedup.clip.bam", "", id)))) %>% # Remove bam suffix.   
    merge(., indvs, by.x = "id", by.y = "fish_ID") %>%             # Add populations.
    mutate(site_full = gsub("[[:space:]]", "", site_full)) %>%     # For consistent merging.
    merge(., sites[,c("Site", "Latitude")], by.x = "site_full",    # Add site latitudes (for plotting).
          by.y = "Site", all.x = T) %>%                            # Keep all individuals.
    filter(site_full %in% shared_pops) %>%                         # Discard unwanted data.
    mutate(site_full = factor(str_sub(gsub("([A-Z])",              # Set site as a factor.
                                      " \\1", site_full), 2)))     # Re-add spaces before second capitals (SanJuan -> San Juan).
  
  df$site_full <- reorder(df$site_full, df$Latitude)               # Order factor wrt latitude.

  return(df)
}

# Obtain and organize PCA information for each dataset.
hdgbs_pruned <- format_eigenvec("../data/pca/hdgbs_pruned.eigenvec")
lcimp_pruned <- format_eigenvec("../data/pca/lcwgs_imputed_pruned.eigenvec")


# Function for visualizing the PCA information for each dataset.
vcf_pca <- \(df, eigenval_file, title, legpos) {

  # Read in eigenvector data and compute %var explained by each PC axis.
  eigvals <- read.delim({{eigenval_file}}, col.names = "eig", header = F)
  PC1 <- sprintf("%0.1f", eigvals[1,1]/sum(eigvals[,1])*100)
  PC2 <- sprintf("%0.1f", eigvals[2,1]/sum(eigvals[,1])*100)
  
  # Assign ggplot object.
  (pc_scatter <- ggplot(df, 
                       aes(x = PC1, y = PC2,
                           group = site_full, 
                           fill = factor(Latitude))) + 
                geom_point(shape = 21, size = 2) + theme_bw() +
                scale_fill_manual(values = c(viridis_pal(option = "D")(length(unique(df$site_full)))),
                                  labels = levels(unique(df$site_full))) +
                theme(legend.title = element_blank(),
                      legend.position = {{legpos}},
                      legend.text = element_text(size = 10)) + 
                ggtitle({{title}}) +
                guides(fill = guide_legend(ncol = 1, byrow = TRUE)) +
                labs(x = paste0("PC1 (", PC1, "%)"), 
                     y = paste0("PC2 (", PC2, "%)"))) 
  
  }

# Visualize PCA for each "vcf-based" dataset. Add legend to first plot only.
(lci_pc <- vcf_pca(df = lcimp_pruned, title = "Imputed lcWGS", legpos = "none",
                   eigenval_file = "../data/pca/lcwgs_imputed_pruned.eigenval"))
(hdg_pr <- vcf_pca(df = hdgbs_pruned, title = "hdGBS", legpos = "right",
                   eigenval_file = "../data/pca/hdgbs_pruned.eigenval"))


# We require a different function for the lcWGS PCA data (it's a *.cov file).
lcwgs_pca <- \(cov_mat, bam_list, title) {
  
  # Read in file/sample names. These must be in the same order as the bam 
  # files used during the PCA calculations!
  # Requires some minor formatting.
  files <- read.table({{bam_list}}, col.names = c("file")) %>% 
    mutate(fish_ID = gsub("309\\-", "309\\.", gsub(".dedup.clip.bam", "", 
                     gsub("^.*IDT_i\\d{1,3}_\\d{1,3}\\.", "", file)))) 
  
  cov    <- read.table({{cov_mat}})               # Read in covariance matrix.
  mmepca <- eigen(cov)                            # Extract eigenvalues.
  eigvec <- as.data.frame(mmepca$vectors) %>%     # Extract eigenvectors.
    `colnames<-`(., gsub("V", "PC", colnames(.))) # Rename columns.

  # % variation explained by each PC axis.
  pca.eigenval.sum <- sum(mmepca$values)
  PC1 <-  sprintf("%0.1f", mmepca$values[1]/pca.eigenval.sum*100)
  PC2 <-  sprintf("%0.1f", mmepca$values[2]/pca.eigenval.sum*100)
  
  pca_dat <- cbind(files, eigvec) %>%            # Bind it all together.
    merge(., indvs, by = "fish_ID") %>%          # Add populations.
    merge(., sites[,c("Site", "Latitude")],      # Add site latitudes.
          by.x = "site_full", 
          by.y = "Site", all.x = T) %>% 
    filter(site_full %in% shared_pops)           # Remove non-shared populations.

    
  # Visualize PCA consistent with above "vcf-based" PCAs.
  (pc_scatter <- ggplot(data  = pca_dat,
                        aes(x = PC1, y = PC2,
                            group = site_full,
                            fill  = factor(Latitude))) +
      geom_point(shape = 21) + theme_bw() +
      scale_fill_manual(values = c(viridis_pal(option = "D")(length(unique(pca_dat$site_full))))) +
      guides(fill = guide_legend(ncol = 1, byrow = TRUE)) +
      theme(legend.title = element_blank(),
            legend.position = "none") +
      labs(x = paste0("PC1 (", PC1, "%)"), 
           y = paste0("PC2 (", PC2, "%)")) +
      ggtitle({{title}}))

}

(lcwgs_full <- lcwgs_pca(cov_mat  = "../data/pca/lcwgs_full_8MSNPs.cov",
                         bam_list = "../data/lcwgs_bam_list_n453.txt",
                         title = "lcWGS") +
    scale_x_continuous(transform = "reverse") +
    scale_y_continuous(transform = "reverse"))


# First, extract the legend of one plot.
pop_legend <- cowplot::get_legend(hdg_pr)

# Arrange PCAs in an organized grid with a shared legend structure.
(all_pcas  <- cowplot::plot_grid(cowplot::plot_grid(plotlist = list(
  hdg_pr + theme(legend.position = "none"), lci_pc, lcwgs_full), 
  ncol = 1, align  = "VH"), pop_legend, 
  rel_widths = c(6,1)))

ggsave("../plots/pcas_ldpruned_main.tiff", dpi = 300,
       width = 8, height = 12, bg = 'white')


################################################################################

scree <- \(eigenval_file, dataset) {
  
  dat <- read_tsv(eigenval_file, col_names = "eigenval") %>% 
    mutate(axis = as.numeric(rownames(.)), PC = paste0("PC", axis)) %>% 
    arrange(axis) %>% mutate(data = dataset)

}

hdg_scr <- scree("../data/pca/hdgbs_full_original.eigenval", dataset = "hdGBS")
imp_scr <- scree("../data/pca/lcwgs_imputed_pruned.eigenval", dataset = "Imputed lcWGS")

full <- rbind(hdg_scr, imp_scr) %>% 
  group_by(data) %>% 
  mutate(var_exp = eigenval/sum(eigenval)) %>% 
  mutate(PC = factor(PC, levels = unique(PC)))

(vars <- ggplot(data = full %>% filter(axis < 10), 
                aes(x = PC, y = var_exp, 
                    colour = data,
                    group = data)) + 
    geom_line(size = 1) + theme_bw() +
    geom_point(shape = 21, colour = "white",
               aes(fill = data), size = 4) +
    labs(x = NULL, y = "Variance explained") +
    scale_y_continuous(labels = scales::percent) +
    theme(legend.title = element_blank(),
          legend.position = c(0.9,0.93),
          legend.box.background = element_rect(color = "black")))

ggsave("../plots/pca_screeplot.tiff", dpi = 300, width = 8, height = 6)

