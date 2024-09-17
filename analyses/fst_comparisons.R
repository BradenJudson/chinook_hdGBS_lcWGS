setwd("~/ots_landscape_genetics/analyses")

library(tidyverse); library(vegan); library(gridExtra); library(pheatmap)


# Ots28 SNP Fsts ---------------------------------------------------------------
# 
# indvs <- read.csv("../landgen_chinook_indvs_corrected.csv", row.names = 1) %>% 
#   mutate(FID = fish_ID) %>% select(c(fish_ID, FID, site_full))
# 
# write.table(indvs, "../data/indv_clusters.txt", 
#             row.names = F, col.names = F, 
#             sep = "\t", quote = F)
# 
# system("plink.exe --vcf ../data/vcfs/hdgbs_subset_134kSNPs_n362_imputed.vcf  --aec --chr NC_056456.1 --fst --within ../data/indv_clusters.txt --out ../data/fst/hdgbsI_Ots28")
# system("plink.exe --vcf ../data/vcfs/hdgbs_subset_134kSNPs_n362_original.vcf --aec --chr NC_056456.1 --fst --within ../data/indv_clusters.txt --out ../data/fst/hdgbsO_Ots28")
# hdgbs_ori <- read.delim("../data/fst/hdgbsO_Ots28.fst") %>% rename("FSToriginal"= FST)
# hdgbs_imp <- read.delim("../data/fst/hdgbsI_Ots28.fst") %>% rename("FSTimputed" = FST)
# 
# hdgbs <- merge(hdgbs_ori[,c(3,5)], hdgbs_imp[,c(3,5)], by = "POS") %>% 
#   filter(FSToriginal > 0)
# 
# 
# ggplot(data = hdgbs,
#        aes(x = FSToriginal,
#            y = FSTimputed)) +
#   geom_point() + theme_bw() +
#   stat_poly_eq(use_label(c("R2", "p")),
#                label.x = "left",
#                label.y = "top",
#                small.p = TRUE)
# 
# 
# 
# 
# 
# 





# Population Fst estimates -----------------------------------------------------

# Retain populations shared by all datasets.
shared_pops <- unique(read.csv("../data/shared_samples_n362.csv")[,"site_full"]) 

# List all pairwise Fst matrix files.
files <- list.files(path = "../data/fst/",
                    pattern = "\\.csv$",
                    full.names = TRUE)

# Read in all pairwise Fst matrices and name them.
fst_mats <- lapply(X = files,
                   FUN = function(x) as.matrix(read.csv(x, 
                                     row.names = 1)) %>% 
                     .[rownames(.) %in% shared_pops,
                       colnames(.) %in% shared_pops]) %>% 
            `names<-`(., str_sub(basename(files), end = -5))

lapply(fst_mats, dim) # Check that dimensions are equal.

# Isolate a consistent clustering order.
pop_order <- hclust(as.dist(fst_mats[[1]]), method = "ward.D2")$order

# Custom heatmap function using pheatmap.
custom_heatmap <- \(dist_mat, title) {
  
  # Set diagonal to NAs.
  diag(dist_mat) <- NA
  
  # Order distance matrix as above and turn off clustering algorithm.
  # Also remove dendrograms and use a light -> dark blue color scheme.
  # Isolate component #4 as that is the plot itself. Helps later.
  pheatmap(mat = dist_mat[pop_order, pop_order] %>% 
             # Add space before capital letter when it follows a lowercase letter.
             `rownames<-`(trimws(str_replace_all(rownames(.), "([A-Z])", " \\1"))) %>% 
             `colnames<-`(trimws(str_replace_all(colnames(.), "([A-Z])", " \\1"))),
           treeheight_row = 0, treeheight_col = 0,
           color = colorRampPalette(c("lightblue1", "dodgerblue3"))(100),
           na_col = NA, legend = F, cluster_rows = F, cluster_cols = F, 
           main = `title`, angle_col = 90)[[4]]
}

# Replace hdGBS imputed with hdGBS non-imputed?

# Make list of relevant plots.
plot_list <- list(
  custom_heatmap(dist_mat = fst_mats[[1]], title = "hdGBS imputed"),
  custom_heatmap(dist_mat = fst_mats[[2]], title = "hdGBS subset"),
  custom_heatmap(dist_mat = fst_mats[[3]], title = "lcWGS imputed")
)

# Save plots over multiple panels. 
heatmap_panels <- do.call(grid.arrange, c(plot_list, ncol = 2, nrow = 2))
ggsave("../plots/fst_heatmap.tiff", heatmap_panels, dpi = 300, width = 15, height = 15)
#### GOAL IS TO MAKE ABOVE 3x2 FOR ALL DATASETS (ncol = 2, nrow = 3) ####


# Compare all pairwise Fst matrices using Mantel statistics.
mantel_list <- list()

for (j in 1:length(fst_mats)) {
  
  for (k in j:length(fst_mats)) {
    
    vm <- vegan::mantel(xdis = fst_mats[[j]],
                        ydis = fst_mats[[k]])
    
    # For organizing outputs.
    df <- data.frame(
      fmat1 = names(fst_mats[j]),
      fmat2 = names(fst_mats[k]),
      mstat = round(vm$statistic, 3),
      msign = vm$signif
    )
    
    mantel_list[[length(mantel_list)+1]] <- df
    
  }
}

# Join output results together in a data frame.
mantel_stats <- do.call("rbind", mantel_list)

# Mantel test matrix for non-subset datasets.
fullfst <- filter(mantel_stats, !grepl("subset", fmat1)) %>% 
  filter(., !grepl("subset", fmat2))

# Mantel test matrix for subset datasets.
subfst  <- filter(mantel_stats,  grepl("subset", fmat1)) %>% 
  filter(.,  grepl("subset", fmat2))

# All relationships are highly significant. 
mean(c(fullfst$msign, subfst$msign)) == 0.001
sd(c(fullfst$msign, subfst$msign)) == 0.000

# Transform the distance matrix into a long-form dataframe.
p2d <- \(df) {
  
  df_mat <- pivot_wider(df[,c(1:3)],
                    names_from = fmat2,
                    values_from = mstat,
                    values_fn = \(x) as.numeric(x[1])) %>% 
    # When values are duplicated, choose the first one.
    column_to_rownames("fmat1") %>% 
    select(names(sort(colSums(is.na(.)), 
                      decreasing = TRUE)))
  
  # Rename and tidy up redundant information.
  dflf <- melt(as.matrix((df_mat))) %>%
    mutate(Var1 = stringr::str_extract(Var1, "[^_]*_[^_]*"),
           Var2 = stringr::str_extract(Var2, "[^_]*_[^_]*"),
           value = na_if(value, 1.0))
}

# Renaming columns essentially 'flips' the diagonal tile geometries where
# the 'full' datasets are above the diagonal and the subset data are below.
ffst <- p2d(fullfst)
sfst <- p2d(subfst) %>% rename("Var2" = Var1, "Var1" = Var2)

# Make a labelling schematic for plotting purposes.
labels <- c("hdGBS", "lcWGS imputed", "lcWGS")

# Plot Mantel statistics between each dataset.
# Shared SNPs and individuals above the diagonal.
# 'Full' datasets below the diagonal.
(mantel_mat <- ggplot() +
  geom_tile(data = ffst,
            aes(Var2, Var1, fill= value)) +
  geom_text(data = ffst, aes(Var2, Var1, label = ifelse(is.na(value), 
            "", sprintf("%0.3f",value))), size = 10, colour = "white") +
  geom_tile(data = sfst, aes(Var2, Var1, fill= value)) +
    geom_text(data = sfst, aes(Var2, Var1, label = ifelse(is.na(value), 
              "", sprintf("%0.3f",value))), size = 10, colour = "white") +
  labs(x = NULL, y = NULL) +
  theme(legend.position = "right",
        panel.grid   = element_blank(),
        legend.ticks = element_blank(),
        legend.text  = element_blank(),
        panel.border = element_blank(),
        axis.ticks   = element_blank(),
        axis.text.x  = element_text(angle = 45, hjust  = 0),
        axis.text = element_text(size = 12, colour = "black")) +
  scale_fill_gradient(low  = "lightblue1", 
                      high = "dodgerblue3", 
                      na.value = NA,
                      name = expression(italic('r'))) +
  scale_x_discrete(position = "top", expand = c(0,0), labels = labels) +
  scale_y_discrete(expand = c(0,0), limits = rev, labels = rev(labels)))
  
ggsave("../plots/fst_mantel_matrix.tiff", dpi = 300,
       width = 10, height = 10)
  

# Need to set up ggsave.
# Make below diagonal common sites and snps.
# Above diagonal is as-is datasets, but among common populations.
# Need to evaluate significance statistics. If all <0.001, state in caption, otherwise
# need to \n and include in each grid cell. 
# Consider better/different color gradient.
# Need improved axis labels. 
# See two links below for using distance matrix in ggplot.
# Need diagonal to be blanks/white.
# For the blue pheatmap plots need to add a space to some names! (i.e., BigQualicum -> Big Qualicum)

# https://stackoverflow.com/questions/48666059/plot-a-re-leveled-pairwise-distance-matrix-in-ggplot2
# https://stackoverflow.com/questions/26838005/putting-x-axis-at-top-of-ggplot2-chart



