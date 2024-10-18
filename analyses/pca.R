setwd("~/ots_landscape_genetics/analyses")

library(tidyverse); library(viridis); library(cowplot)

indvs <- read.csv("../data/landgen_chinook_indvs.csv", na.strings = "") %>% 
  mutate(site_full = gsub("[[:space:]]Brood", "", site_full))
indvs[indvs$site_full == "Upper Pitt", "site_full"] <- "Pitt"

# Retain populations shared by all datasets.
shared_pops <- unique(read.csv("../data/shared_samples_n362.csv")[,"site_full"]) 
sites <- read.delim(file = "../data/ch2023_sequenced.txt") %>% 
  mutate(Site = gsub("[[:space:]]", "", Site),
         sitenum = as.numeric(rownames(.))) %>% 
  arrange(Latitude) %>% 
  filter(Site %in% shared_pops)

# --pca 1000 just ensures that all possible PCs are calculated (so %var explained is accurate).
system("plink.exe --vcf ../data/vcfs/hdgbs_full_maf5_original.vcf.gz --double-id --aec --pca 1000 --out ../data/pca/hdgbs_full_original")
system("plink.exe --vcf ../data/vcfs/lcWGS_full_8MSNPs_imputed.vcf.gz --double-id --aec --pca 1000 --out ../data/pca/lcWGS_full_8MSNPs_imputed")
system("plink.exe --vcf ../data/vcfs/hdgbs_subset_134kSNPs_n362_original.vcf.gz --double-id --aec --pca 1000 --out ../data/pca/hdgbs_subset_original")
system("plink.exe --vcf ../data/vcfs/lcwgs_subset_134kSNPs_imputed.vcf.gz --double-id --aec --pca 1000 --out ../data/pca/lcWGS_subset_imputed")

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
hdgbs_full_pca  <- format_eigenvec("../data/pca/hdgbs_full_original.eigenvec")
hdgbs_sub_pca   <- format_eigenvec("../data/pca/hdgbs_subset_original.eigenvec")
lcwgs_s.imputed <- format_eigenvec("../data/pca/lcWGS_full_8MSNPs_imputed.eigenvec")
lcwgs_f.imputed <- format_eigenvec("../data/pca/lcWGS_subset_imputed.eigenvec")

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
                geom_point(shape = 21) + theme_bw() +
                scale_fill_manual(values = c(viridis_pal(option = "D")(length(unique(df$site_full)))),
                                  labels = levels(unique(df$site_full))) +
                theme(legend.title = element_blank(),
                      legend.position = {{legpos}},
                      legend.text = element_text(size = 10)) + 
                ggtitle({{title}}) +
                guides(fill = guide_legend(ncol = 2, byrow = TRUE)) +
                labs(x = paste0("PC1 (", PC1, "%)"), 
                     y = paste0("PC2 (", PC2, "%)"))) 
  
  }

# Visualize PCA for each "vcf-based" dataset. Add legend to first plot only.
(hdg_f <- vcf_pca(df = hdgbs_full_pca, "../data/pca/hdgbs_full_original.eigenval", title = "hdGBS", legpos = "right") +
    scale_x_continuous(transform = "reverse") + scale_y_continuous(transform = "reverse"))
(hdg_s <- vcf_pca(df = hdgbs_sub_pca,  "../data/pca/hdgbs_subset_original.eigenval", title = "hdGBS subset", legpos = "none") +
    scale_x_continuous(transform = "reverse"))
(lci_f <- vcf_pca(df = lcwgs_f.imputed, "../data/pca/lcWGS_full_8MSNPs_imputed.eigenval", title = "lcWGS imputed", legpos = "none"))
(lci_s <- vcf_pca(df = lcwgs_s.imputed, "../data/pca/lcWGS_subset_imputed.eigenval", title = "lcWGS imputed subset", legpos = "none"))

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
      guides(fill = guide_legend(ncol = 2, byrow = TRUE)) +
      theme(legend.title = element_blank(),
            legend.position = "none") +
      labs(x = paste0("PC1 (", PC1, "%)"), 
           y = paste0("PC2 (", PC2, "%)")) +
      ggtitle({{title}}))

}

(lcwgs_full <- lcwgs_pca(cov_mat  = "../data/pca/lcwgs_full_8MSNPs.cov",
                         bam_list = "../data/lcwgs_bam_list_n453.txt",
                         title = "lcWGS"))

(lcwgs_subs <- lcwgs_pca(cov_mat  = "../data/pca/angsd_subset134ksnps.cov",
                         bam_list = "../data/lcwgs_bam_list_n453.txt",
                         title = "lcWGS subset") + 
                         scale_x_continuous(transform = "reverse"))

# First, extract the legend of one plot.
pop_legend <- cowplot::get_legend(hdg_f)

# Arrange PCAs in an organized grid with a shared legend structure.
(all_pcas  <- cowplot::plot_grid(cowplot::plot_grid(plotlist = list(
  hdg_f + theme(legend.position = "none"), hdg_s, lci_f, lci_s, 
  lcwgs_full, lcwgs_subs), ncol = 2, align  = "VH"), pop_legend, 
  ncol = 2, rel_widths = c(6,1)))

ggsave("../plots/all_pcas.tiff", dpi = 300,
       width = 12, height = 12, bg = 'white')
