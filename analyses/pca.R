setwd("~/ots_landscape_genetics/analyses")

library(tidyverse); library(viridis); library(cowplot); library(ggpubr)

indvs <- read.csv("../data/landgen_chinook_indvs.csv", na.strings = "") %>% 
  mutate(site_full = gsub("[[:space:]]Brood", "", site_full))
indvs[indvs$site_full == "Upper Pitt", "site_full"] <- "Pitt"

# Retain populations shared by all datasets.
shared_pops <- unique(read.csv("../data/shared_samples_n361.csv")[,"site_full"]) 
sites <- readxl::read_xlsx("../data/chinook_lineages.xlsx", sheet = 2) %>% 
  mutate(Site = gsub("[[:space:]]", "", Site)) %>% 
  arrange(Latitude) %>% 
  mutate(sitenum = as.numeric(rownames(.)))


# --pca 1000 just ensures that all possible PCs are calculated (so %var explained is accurate).
system("plink.exe --vcf ../data/vcfs_n385/hdgbs_maf5_n385_ldpruned.vcf.gz  --aec --double-id --pca 1000 --out ../data/pca_n385/hdgbs_pruned")
system("plink.exe --vcf ../data/vcfs_n385/chinook_lcwgs_maf005_n385_imputed_ldpruned.vcf.gz --aec --double-id --pca 1000 --out ../data/pca_n385/lcwgs_imputed_pruned")
# system("plink.exe --vcf ../data/vcfs/hdgbs_indvsFiltered.vcf.gz --aec --double-id --pca 1000 --out ../data/pca/hdgbs_indvsFiltered")


# Format PCA output data from above.
format_eigenvec <- \(eigenvec_file) {
  df <- read.table(eigenvec_file) %>% .[,2:ncol(.)] %>%            # Remove unnecessary first column.
    `colnames<-`(., c("id", paste0("PC", 1:(ncol(.)-1)))) %>%      # Rename columns.
    mutate(id = gsub("309\\-", "309\\.", gsub('^.*\\_[0-9.]+', "", # Make label formatting consistent.
                gsub(pattern = "\\.dedup.clip.bam", "", id)))) %>% # Remove bam suffix.   
    merge(., indvs, by.x = "id", by.y = "fish_ID") %>%             # Add populations.
    mutate(site_full = gsub("[[:space:]]", "", site_full)) %>%     # For consistent merging.
    merge(., sites[,c("Site", "Latitude", "Region")], 
          by.x = "site_full",                                      # Add site latitudes (for plotting).
          by.y = "Site", all.x = T) %>%                            # Keep all individuals.
    # filter(site_full %in% shared_pops) %>%                       # Discard unwanted data.
    mutate(Region = factor(Region),
           site_full = factor(str_sub(gsub("([A-Z])",              # Set site as a factor.
                                      " \\1", site_full), 2)))     # Re-add spaces before second capitals (SanJuan -> San Juan).
  
  df$Region <- reorder(df$Region, df$Latitude)                     # Order factor wrt latitude.

  return(df)
}

# Obtain and organize PCA information for each dataset.
hdgbs_pruned <- format_eigenvec("../data/pca_n385/hdgbs_pruned.eigenvec")
lcimp_pruned <- format_eigenvec("../data/pca_n385/lcwgs_imputed_pruned.eigenvec")

# Summarize variation within populations. 
hd_sd <- hdgbs_pruned %>% group_by(site_full) %>% 
  summarise(pc1sd = sd(PC1), pc2sd = sd(PC2)) %>% 
  mutate(dataset = "hdGBS")
lc_sd <- lcimp_pruned %>% group_by(site_full) %>% 
  summarise(pc1sd = sd(PC1), pc2sd = sd(PC2)) %>% 
  mutate(dataset = "Imputed lcWGS")

sds <- rbind(hd_sd, lc_sd) %>% 
  pivot_longer(cols = c(pc1sd, pc2sd)) %>% 
  mutate(PC = case_when(name == "pc1sd" ~ "PC1",
                        name == "pc2sd" ~ "PC2"))

ggboxplot(data = sds, x = "dataset", 
          y = "value", facet.by = "PC") +
  stat_compare_means(comparisons = list(c("hdGBS", "Imputed lcWGS")),
                     method = "t.test", paired = T) +
  labs(x = NULL, y = "Standard deviation") +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA, colour = NA))


ggsave("../plots/pca_stdevs.tiff", dpi = 300,
       width = 10, height = 5)

t.test(hd_sd$pc1sd, lc_sd$pc1sd, paired = T)
mean(hd_sd$pc1sd); mean(lc_sd$pc1sd)
t.test(hd_sd$pc2sd, lc_sd$pc2sd, paired = T)
mean(hd_sd$pc2sd); mean(lc_sd$pc2sd)


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
                           fill = factor(Region),
                           shape = factor(Region))) + 
                geom_point(size = 2) + theme_bw() +
                scale_shape_manual(values = c(rep(c(21,23),8))) +
                theme(legend.title = element_blank(),
                      legend.position = {{legpos}},
                      legend.text = element_text(size = 10)) + 
                ggtitle({{title}}) +
                guides(shape = guide_legend(nrow = 3, byrow = TRUE)) +
                labs(x = paste0("PC1 (", PC1, "%)"), 
                     y = paste0("PC2 (", PC2, "%)"))) 
  
  }

# Visualize PCA for each "vcf-based" dataset. Add legend to first plot only.
(lci_pc <- vcf_pca(df = lcimp_pruned, title = "Imputed lcWGS", legpos = "none",
                   eigenval_file = "../data/pca_n385/lcwgs_imputed_pruned.eigenval") +
    scale_y_continuous(transform = "reverse"))
(hdg_pr <- vcf_pca(df = hdgbs_pruned, title = "hdGBS", legpos = "right",
                   eigenval_file = "../data/pca/hdgbs_pruned.eigenval") +
    scale_x_continuous(transform = "reverse"))




# (SG <- rbind(hdg_pr$data %>% mutate(dataset = "hdGBS"),
#            lci_pc$data %>% mutate(dataset = "Imputed lcWGS")) %>% 
#   dplyr::rename(Population = site_full) %>% 
#   filter(Region %in% c("Vancouver Island (east coast)")) %>%
#   ggplot(data = ., aes(x = PC1, y = PC2,
#                        color = Population)) +
#   geom_polygon(stat = "ellipse", aes(fill = Population), alpha = 1/5) +
#   geom_point(aes(shape = Population)) +
#   facet_wrap(~dataset, scales = "free") +
#   scale_shape_manual(values = seq(1:7)) +
#   theme_bw() +
#   theme(legend.position = "top",
#         legend.title = element_blank()))
# 
# 
# reg_list <- list()
# 
# for (region in c("Stikine/Taku Rivers", "Yukon River", "Fraser River", "Thompson River",
#                  "Vancouver Island (west coast)", "Vancouver Island (east coast)")) {
#   
#   reg_list[[region]] <- rbind(hdg_pr$data %>% mutate(dataset = "hdGBS"),
#                               lci_pc$data %>% mutate(dataset = "Imputed lcWGS")) %>% 
#     dplyr::rename(Population = site_full) %>% 
#     filter(Region %in% region) %>%
#     ggplot(data = ., aes(x = PC1, y = PC2,
#                          color = Population)) +
#     geom_polygon(stat = "ellipse", aes(fill = Population), alpha = 1/5) +
#     geom_point(aes(shape = Population)) +
#     facet_wrap(~dataset, scales = "free") +
#     scale_shape_manual(values = seq(1:11)) +
#     theme_bw() +
#     theme(legend.position = "top",
#           legend.title = element_blank())
#   
# }
# 
# (reg_pca <- cowplot::plot_grid(plotlist = reg_list))
# ggsave("../plots/supp_figs/regional_pcas.tiff", dpi = 300, height = 10, width = 15)



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
    merge(., sites[,c("Site", "Latitude",        # Add site latitudes. 
                      "Region")],       
          by.x = "site_full", 
          by.y = "Site", all.x = T) 

  pca_dat$Region <- reorder(pca_dat$Region, pca_dat$Latitude)
    
  # Visualize PCA consistent with above "vcf-based" PCAs.
  (pc_scatter <- ggplot(data  = pca_dat,
                        aes(x = PC1, y = PC2,
                            group = site_full,
                            fill  = factor(Region),
                            shape = factor(Region))) +
      geom_point(size = 2) + theme_bw() +
      scale_shape_manual(values = c(rep(c(21,23),8))) +
      # scale_fill_manual(values = c(viridis_pal(option = "D")(length(unique(pca_dat$site_full))))) +
      guides(fill = guide_legend(ncol = 1, byrow = TRUE)) +
      theme(legend.title = element_blank(),
            legend.position = "none") +
      labs(x = paste0("PC1 (", PC1, "%)"), 
           y = paste0("PC2 (", PC2, "%)")) +
      ggtitle({{title}}))

}

(lcwgs_full <- lcwgs_pca(cov_mat  = "../data/pca_n385/angsd_n385_ldpruned.cov",
                         bam_list = "../data/bam_list_n385.txt",
                         title = "lcWGS") +
    scale_x_continuous(transform = "reverse"))


# First, extract the legend of one plot.
pop_legend <- cowplot::get_legend(hdg_pr)


# Arrange PCAs in an organized grid with a shared legend structure.
(all_pcas  <- cowplot::plot_grid(cowplot::plot_grid(plotlist = list(
  hdg_pr + theme(legend.position = "none"), lci_pc, lcwgs_full), 
  ncol = 3, align  = "VH"), pop_legend, ncol = 1, rel_heights = c(5,1)))



ggsave("../plots/pcas_ldpruned_main2.tiff", dpi = 300,
       width = 14, height = 5.5, bg = 'white')


################################################################################

scree <- \(eigenval_file, dataset) {
  
  dat <- read_tsv(eigenval_file, col_names = "eigenval") %>% 
    mutate(axis = as.numeric(rownames(.)), PC = paste0("PC", axis)) %>% 
    arrange(axis) %>% mutate(data = dataset)

}

hdg_scr <- scree("../data/pca/hdgbs_full_original.eigenval", dataset = "hdGBS")
imp_scr <- scree("../data/pca/lcwgs_imputed_pruned.eigenval", dataset = "Imputed lcWGS")
lcw_scr <- data.frame(eigenval = eigen(read.table("../data/pca_n385/angsd_n385_ldpruned.cov"))$values) %>% 
  mutate(axis = as.numeric(rownames(.)), PC = paste0("PC", axis)) %>% 
  arrange(axis) %>% mutate(data = "lcWGS")

full <- rbind(hdg_scr, imp_scr, lcw_scr) %>% 
  group_by(data) %>% 
  mutate(var_exp = eigenval/sum(eigenval)) %>% 
  mutate(PC = factor(PC, levels = unique(PC)))

(vars <- ggplot(data = full %>% filter(axis < 10), 
                aes(x = PC, y = var_exp, 
                    colour = data,
                    group = data)) + 
    geom_line(linewidth = 1) + theme_bw() +
    geom_point(shape = 21, colour = "white",
               aes(fill = data), size = 4) +
    labs(x = NULL, y = "Variance explained") +
    scale_y_continuous(labels = scales::percent) +
    theme(legend.title = element_blank(),
          legend.position = c(0.9,0.91),
          legend.box.background = element_rect(color = "black")) +
    annotate("segment", y = 0.04, yend = .015, x = "PC4", xend = "PC4",
             arrow = arrow(type = "closed", length = unit(0.02, "npc"))))

ggsave("../plots/supp_figs/pca_screeplot.tiff", dpi = 300, width = 8, height = 6)

