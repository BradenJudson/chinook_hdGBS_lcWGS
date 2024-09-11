setwd("~/ots_landscape_genetics/analyses")

library(tidyverse); library(vegan); library(ggpmisc); library(pheatmap)


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

shared_pops <- unique(read.csv("../data/shared_samples_n362.csv")[,"site_full"]) 

files <- list.files(path = "../data/fst/",
                    pattern = "\\.csv$",
                    full.names = TRUE)

fst_mats <- lapply(X = files,
                   FUN = function(x) as.matrix(read.csv(x, 
                                     row.names = 1)) %>% 
                     .[rownames(.) %in% shared_pops,
                       colnames(.) %in% shared_pops]) %>% 
            `names<-`(., str_sub(basename(files), end = -5))

lapply(fst_mats, dim)

# f1 <- melt(t(fst_mats[[2]]), value.name = "f1")
# f2 <- melt(as.matrix(fst_mats[[1]]), value.name = "f2")
(j <- vegan::mantel(xdis = fst_mats[[1]], ydis = fst_mats[[2]]))
# 
# j <- merge(f1, f2, by = c("Var1", "Var2")) %>% 
#   mutate(f1 = f1*rnorm(n = nrow(.), sd = 1/6, mean = 1))
# 
# ggplot(data = j,
#        aes(x = f1, y = f2)) +
#   geom_abline(slope = 1) +
#   geom_smooth(method = "lm", se = F) +
#   geom_point(alpha = 1/2) +
#   theme_bw() +
#   theme(panel.grid = element_blank()) +
#   stat_poly_eq(use_label(c("R2", "p")),
#                label.x = "left",
#                label.y = "top",
#                small.p = TRUE) +
#   labs(x = expression(F[ST]),
#        y = expression(F[ST]))

# Isolate a consistent clustering order.
pop_order <- hclust(as.dist(fst_mats[[1]]))$order

# Custom heatmap function using pheatmap.
custom_heatmap <- \(dist_mat) {
  
  # Set diagonal to NAs.
  diag(dist_mat) <- NA
  
  # Order distance matrix as above and turn off clustering algorithm.
  # Also remove dendrograms and use a light -> dark blue color scheme.
  # Isolate component #4 as that is the plot itself. Helps later.
  pheatmap(mat = dist_mat[pop_order, pop_order],
           treeheight_row = 0, treeheight_col = 0,
           color = colorRampPalette(c("lightblue1", "dodgerblue3"))(100),
           na_col = NA, legend = F, cluster_rows = F, cluster_cols = F)[[4]]
}

# make list of relevant plots.
plot_list <- list(
  custom_heatmap(dist_mat = fst_mats$hdgbs_subset_134kSNPs_n362_original),
  custom_heatmap(dist_mat = fst_mats$hdgbs_snps_maf005_imputed)
)

# Save plots over multiple panels. 
heatmap_panels <- do.call(grid.arrange, plot_list)
ggsave("../plots/fst_heatmap.tiff", heatmap_panels, dpi = 300, width = 10, height = 16)
