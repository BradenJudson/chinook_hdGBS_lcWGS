setwd("~/ots_landscape_genetics/analyses")

library(tidyverse); library(vegan); library(gridExtra); library(pheatmap)
library(reshape2); library(data.table); library(ggpmisc)


# Site-wise Fst Estimates ------------------------------------------------------

# Identify, read in and reformat site-wise Fst files.
snp_fst <- list.files("../data/fst/", pattern = "weir.fst", full.names = T)
fst_dat <- lapply(snp_fst,
                  FUN = \(x) fread(x) %>% 
                        rename("Fst" = 3)) %>% 
  # Rename list elements and cut off string ".weir.fst".
  `names<-`(., str_sub(basename(snp_fst), end = -10))

# Function to collapse list into a wide format dataframe with better colnames.
wfsnps <- \(x) map_df(x, ~as.data.frame(.x), .id = "dataset") %>% 
  rename("Fst" = 4) %>% pivot_wider(., names_from = dataset, values_from = Fst)

# Return data for subset and "full" datasets considered separately.
sub_fst_snps  <- wfsnps(fst_dat[grepl("subset",  names(fst_dat))])
full_fst_snps <- wfsnps(fst_dat[!grepl("subset", names(fst_dat))])

# Establish improved labelling vector.
plot_labs <- c("chinook_imputed_8M" = "Imputed lcWGS",
               "hdgbs_full_imputed" = "hdGBS imputed",
               "hdgbs_subset_134kSNPs_imputed" = "hdGBS subset imputed",
               "lcwgs_imputed_134kSNPs_subset" = "lcWGS subset imputed")

# From the wide-form dataframe, create a scatter plot of 
# estimated site-wise Fst values. Set constant axis boundaries/labels.
# Also print R2 and show OLS best fit line.
scatterFST <- function(df, x_axis, y_axis) {
  ggplot(data  = df,
         aes(x = {{x_axis}},
             y = {{y_axis}})) +
    theme_classic() +
    geom_point(shape = 21,
               fill  = "gray",
               colour= "black",
               alpha = 3/4) +
    stat_smooth(method = "lm",
                alpha  = 1/6,
                colour = "black",
                linewidth = 2,
                lineend = "round") +
    stat_smooth(method = "lm",
                alpha  = 1/6,
                color  = "gray90",
                lineend = "round") +
    scale_x_continuous(limits = c(-0.05,1),
                       breaks = seq(-0.1, 1, 1/4)) +
    scale_y_continuous(limits = c(-0.05,1),
                       breaks = seq(-0.1, 1, 1/4)) +
    stat_poly_eq(use_label(c("R2", "p")),
                 label.x = "left",
                 label.y = "top",
                 small.p = TRUE) +
    labs(x = plot_labs[[deparse(substitute(x_axis))]],
         y = plot_labs[[deparse(substitute(y_axis))]])
}

# Subset data scatter plots.
(s1 <- scatterFST(sub_fst_snps, hdgbs_subset_134kSNPs_imputed, lcwgs_imputed_134kSNPs_subset))

# Full data scatter plots.
(f1 <- scatterFST(full_fst_snps, hdgbs_full_imputed, chinook_imputed_8M))

# Function for creating Manhattan plots for Ots28 only. 
manhat3 <- \(df) {
  
  # First pivot from wide to long form and choose chromosome.
  lf_df <- pivot_longer(data = df[df$CHROM == "NC_056456.1",], 
                        cols = 3:ncol(df),
                        values_to = "Fst",
                        names_to  = "dataset")
  
  ggplot(data  = lf_df, 
         aes(x = POS/1e6, 
             y = Fst)) +
    geom_point(shape = 21,
               fill  = "gray80",
               alpha = 4/5) +
    labs(x = "Position (Mbp)",
         y = expression(F[ST])) +
    theme_bw() +
    facet_wrap( ~ dataset, ncol = 1, scales = "free_x", 
                labeller = as_labeller(plot_labs)) +
    theme(strip.background = element_rect(color = NA, fill = NA),
          plot.background  = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          strip.placement = "inside",
          strip.text.x = element_text(size = 12),
          axis.line    = element_line(colour = "black")) 
    # free_x necessary to keep bottom panel across facets.
    
}

# Arrange Manhattan and scatter plots together.
# Left (A) is Ots28 only, and right (B) is genome-wide.
cowplot::plot_grid(
  manhat3(df = sub_fst_snps),
  cowplot::plot_grid(s1, f1, 
           ncol = 1, scale = 9/10),
  labels = c("A)", "B)"),
  rel_widths = c(2, 1)
)

ggsave("../plots/site_wise_fst.tiff", 
       dpi = 300, width = 14, height = 7)
ggsave("../plots/site_wise_fst.png", 
       dpi = 300, width = 14, height = 7)


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
  # Also remove dendrograms and use a light -> dark blue colour scheme.
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

# Make list of relevant plots.
plot_list <- list(
  custom_heatmap(dist_mat = fst_mats$hdgbs_original_full, title = "hdGBS"),
  custom_heatmap(dist_mat = fst_mats$hdgbs_original_134kSNPs_n362_subset, title = "hdGBS subset"),
  custom_heatmap(dist_mat = fst_mats$lcwgs_imputed_8MSNPs_full, title = "lcWGS imputed"),
  custom_heatmap(dist_mat = fst_mats$lcwgs_imputed_134kSNPs_subset, title = "lcWGS imputed subset"),
  custom_heatmap(dist_mat = fst_mats$lcwgs_angsd_weightedfst_matrix, title = "lcWGS"),
  custom_heatmap(dist_mat = fst_mats$lcwgs_angsd_subset134kSNPs_fst_matrix, title = "lcWGS subset")
)

# Save plots over multiple panels. 
heatmap_panels <- do.call(grid.arrange, c(plot_list, ncol = 2))
ggsave("../plots/fst_heatmap.tiff", heatmap_panels, dpi = 300, width = 15, height = 20)
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
    
    mantel_list[[length(mantel_list)+1]] <- df }
  
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
labels <- c("hdGBS", "lcWGS","lcWGS imputed")

# Plot Mantel statistics between each dataset.
# Shared SNPs and individuals above the diagonal.
# 'Full' datasets below the diagonal.
(mantel_mat <- ggplot() +
  geom_tile(data = ffst,
            aes(Var2, Var1, fill= value)) +
  geom_text(data = ffst, aes(Var2, Var1, label = ifelse(is.na(value), 
            "", sprintf("%0.3f",value))), size = 10, colour = "white") +
  geom_tile(data = sfst, 
            aes(Var2, Var1, fill= value)) +
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
        axis.text    = element_text(size = 12, colour = "black")) +
  scale_fill_gradient(low  = "lightblue1", 
                      high = "dodgerblue3", 
                      na.value = NA,
                      name = expression(italic('  r'))) +
  scale_x_discrete(position = "top", expand = c(0,0), labels = labels) +
  scale_y_discrete(expand = c(0,0), labels = rev(labels), limits = rev))
  
ggsave("../plots/fst_mantel_matrix.tiff", dpi = 300,
       width = 10, height = 10)

