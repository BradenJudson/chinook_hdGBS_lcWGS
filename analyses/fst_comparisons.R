setwd("~/ots_landscape_genetics/analyses")

library(tidyverse); library(vegan); library(gridExtra)
library(ggpmisc); library(cowplot); library(pheatmap)


# Population Fst estimates -----------------------------------------------------

hd_fst <- read.delim("../data/fst_n385/hdgbs_n385_hudson_fst.fst.summary", sep = "") 
im_fst <- read.delim("../data/fst_n385/lcwgs_n385_imputed_pruned_fst.fst.summary", sep = "") 

lc_fst <- read.delim("../data/fst_n385/lcwgs_fst_estimates_n385.txt", sep = "") %>% 
  mutate(pops = gsub("san_","san",gsub("big_","big", 
                gsub("upper_", "upper", gsub(".globalFst","",
                gsub("14_fst/global_fst_n385/", "", file))))),
         POP1 = gsub("Upperpitt","UpperPitt", gsub("Bigqualicum","BigQualicum", 
                gsub("Sanjuan","SanJuan", gsub("Bigsalmon","BigSalmon",
                tools::toTitleCase(gsub("_[^_]+$", "", pops)))))),
         POP2 = gsub("Upperpitt","UpperPitt", gsub("Bigqualicum","BigQualicum", 
                gsub("Sanjuan","SanJuan", gsub("Bigsalmon","BigSalmon",
                tools::toTitleCase(gsub(".*\\_", "", pops))))))) %>% 
  select(c("POP1", "POP2", "FstWeighted"))

(fst_lm1 <- lm(im_fst$HUDSON_FST ~ hd_fst$HUDSON_FST)); summary(fst_lm1)
(fst_lm2 <- lm(lc_fst$FstWeighted ~ hd_fst$HUDSON_FST)); summary(fst_lm2)

fst <- merge(hd_fst, im_fst, by = c("X.POP1", "POP2")) %>% 
  rename(POP1 = X.POP1) %>% 
  merge(., lc_fst, c("POP1", "POP2")) %>% 
  rename(hdgbs = HUDSON_FST.x, lcimp = HUDSON_FST.y,
         lcwgs = FstWeighted)

t.test(fst$hdgbs, fst$lcimp, paired = T)
t.test(fst$hdgbs, fst$lcwgs, paired = T)
t.test(fst$lcimp, fst$lcwgs, paired = T)

fct_labs <- c("lcwgs" = "lcWGS", "lcimp" = "Imputed lcWGS")

(pop_fst <- fst %>% 
    pivot_longer(cols = c("lcimp", "lcwgs")) %>% 
   ggplot(data = ., aes(x = hdgbs, y = value)) + 
   geom_point(shape = 21, fill = "gray85", size = 1.5, alpha = 2/3) +
   theme_bw() + 
   labs(x = expression(`hdGBS  F`[ST]), 
        y = expression(F[ST])) +
   geom_abline(slope = 1, intercept = 0, colour = "blue2", linewidth = 1) +
   stat_smooth(method = "lm", se = F, colour = "black") +
    stat_poly_eq(use_label(c("R2")), label.x = "right",
                 label.y = "bottom") +
    facet_wrap(~name,labeller = as_labeller(fct_labs)) +
    theme(strip.background = element_blank()))

ggsave("../plots/supp_figs/fst_lm.tiff", dpi = 300, width = 10, height = 6)

# Convert longform distances into pairwise distance matrices.
hfstmat <- ConGenFunctions::df2pmat(df = hd_fst, var = "HUDSON_FST",  grp1 = "X.POP1", grp2 = "POP2")
ifstmat <- ConGenFunctions::df2pmat(df = im_fst, var = "HUDSON_FST",  grp1 = "X.POP1", grp2 = "POP2")
lfstmat <- ConGenFunctions::df2pmat(df = lc_fst, var = "FstWeighted", grp1 = "POP1",   grp2 = "POP2")

fst_mats <- setNames(list(hfstmat, ifstmat, lfstmat), 
                     c("hdGBS", "imputedlcWGS", "lcWGS"))

lapply(fst_mats, dim) # Check that dimensions are equal.

# Isolate a consistent clustering order.
pop_order <- hclust(as.dist(fst_mats[[1]]), method = "ward.D2")$order

# Custom heatmap function using pheatmap.
custom_heatmap <- \(dist_mat, title, LR) {
  
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
           main = `title`, angle_col = 90, labels_row = LR)[[4]]
  
}


plot_list <- list(
  custom_heatmap(dist_mat = fst_mats$hdGBS, title = "hdGBS", LR = ""),
  custom_heatmap(dist_mat = fst_mats$imputedlcWGS, title = "Imputed lcWGS", LR = ""),
  custom_heatmap(dist_mat = fst_mats$lcWGS, title = "lcWGS", LR = NULL)
)

# Save plots over multiple panels. 
(heatmap_panels <- plot_grid(plotlist = plot_list, nrow = 1, rel_widths = c(1, 1, 1.12)))
ggsave("../plots/fst_heatmap_n385.tiff", heatmap_panels, dpi = 300, width = 22, height = 8, bg = 'white')

# Compare all pairwise Fst matrices using Mantel statistics.
# # mantel_list <- list()
# 
# for (j in 1:length(fst_mats)) {
# 
#   for (k in j:length(fst_mats)) {
# 
#     vm <- vegan::mantel(xdis = fst_mats[[j]],
#                         ydis = fst_mats[[k]])
# 
#     # For organizing outputs.
#     df <- data.frame(
#       fmat1 = names(fst_mats[j]),
#       fmat2 = names(fst_mats[k]),
#       mstat = round(vm$statistic, 3),
#       msign = vm$signif
#       )
# 
#     mantel_list[[length(mantel_list)+1]] <- df }
# 
# }
# 
# # Join output results together in a data frame.
# mantel_stats <- do.call("rbind", mantel_list)
# 
# # Mantel test matrix for non-subset datasets.
# fullfst <- filter(mantel_stats, !grepl("subset", fmat1)) %>%
#   filter(., !grepl("subset", fmat2))
# 
# # Mantel test matrix for subset datasets.
# subfst  <- filter(mantel_stats,  grepl("subset", fmat1)) %>%
#   filter(.,  grepl("subset", fmat2))
# 
# # All relationships are highly significant.
# mean(c(fullfst$msign, subfst$msign)) == 0.001
# sd(c(fullfst$msign, subfst$msign)) == 0.000
# 
# # Transform the distance matrix into a long-form dataframe.
# p2d <- \(df) {
# 
#   df_mat <- pivot_wider(df[,c(1:3)],
#                     names_from = fmat2,
#                     values_from = mstat,
#                     values_fn = \(x) as.numeric(x[1])) %>%
#     # When values are duplicated, choose the first one.
#     column_to_rownames("fmat1") %>%
#     select(names(sort(colSums(is.na(.)),
#                       decreasing = TRUE)))
# 
#   # Rename and tidy up redundant information.
#   dflf <- melt(as.matrix((df_mat))) %>%
#     mutate(Var1 = stringr::str_extract(Var1, "[^_]*_[^_]*"),
#            Var2 = stringr::str_extract(Var2, "[^_]*_[^_]*"),
#            value = na_if(value, 1.0))
# }
# 
# # Renaming columns essentially 'flips' the diagonal tile geometries where
# # the 'full' datasets are above the diagonal and the subset data are below.
# ffst <- p2d(fullfst)
# sfst <- p2d(subfst) %>% rename("Var2" = Var1, "Var1" = Var2)
# 
# # Make a labelling schematic for plotting purposes.
# labels <- c("hdGBS", "lcWGS","lcWGS imputed")
# 
# # Plot Mantel statistics between each dataset.
# # Shared SNPs and individuals above the diagonal.
# # 'Full' datasets below the diagonal.
# (mantel_mat <- ggplot() +
#   geom_tile(data = ffst,
#             aes(Var2, Var1, fill= value)) +
#   geom_text(data = ffst, aes(Var2, Var1, label = ifelse(is.na(value),
#             "", sprintf("%0.3f",value))), size = 10, colour = "white") +
#   geom_tile(data = sfst,
#             aes(Var2, Var1, fill= value)) +
#   geom_text(data = sfst, aes(Var2, Var1, label = ifelse(is.na(value),
#             "", sprintf("%0.3f",value))), size = 10, colour = "white") +
#   labs(x = NULL, y = NULL) +
#   theme(legend.position = "right",
#         panel.grid   = element_blank(),
#         legend.ticks = element_blank(),
#         legend.text  = element_blank(),
#         panel.border = element_blank(),
#         axis.ticks   = element_blank(),
#         axis.text.x  = element_text(angle = 45, hjust  = 0),
#         axis.text    = element_text(size = 12, colour = "black")) +
#   scale_fill_gradient(low  = "lightblue1",
#                       high = "dodgerblue3",
#                       na.value = NA,
#                       name = expression(italic('  r'))) +
#   scale_x_discrete(position = "top", expand = c(0,0), labels = labels) +
#   scale_y_discrete(expand = c(0,0), labels = rev(labels), limits = rev))
# 
# ggsave("../plots/fst_mantel_matrix.tiff", dpi = 300,
#        width = 10, height = 10)


