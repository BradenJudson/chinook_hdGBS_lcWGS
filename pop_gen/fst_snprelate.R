setwd("~/ots_landscape_genetics/pop_gen")

library(tidyverse); library(SNPRelate); library(svMisc)

# Close GDS envirnoment (if not starting fresh, this can cause issues with writing).
gdsfmt::showfile.gds(closeall = T)

# Convert VCF to GDS object.
snpgdsVCF2GDS(vcf.fn = "../data/vcfs/hdgbs_snps_maf005_imputed.vcf.gz", out.fn = "hdgbs_gds.gds")

# Read in GDS object.
(genofile <- snpgdsOpen("hdgbs_gds.gds"))

# Read in individual sample information.
indvs <- read.csv("../data/landgen_chinook_indvs_corrected.csv", row.names = 1)

# Create a dataframe of all possible combinations of populations.
pop_pairs <- t(as.data.frame(combn(unique(indvs$site_full), 2)))

# List for populating.
flist <- list()

# Calculates FST between each pair of populations.
for (i in 1:nrow(pop_pairs)) {
  
  pops <- pop_pairs[i,]
  
  samples <- indvs[indvs$site_full %in% c(pops), "fish_ID"]
  pop_vec <- indvs[indvs$site_full %in% c(pops), "site_full"]
  
  gfst <- snpgdsFst(genofile, 
                 sample.id  = c(samples),
                 # Population must be a factor!
                 population = as.factor(pop_vec),
                 autosome.only = FALSE,
                 method  = "W&C84",
                 verbose = FALSE)
  
  # Arrange outputs in a small dataframe.
  fst_dat <- data.frame(
    pop1 = pops[1],
    pop2 = pops[2],
    meanFST = round(gfst$MeanFst, 5),
    weightedFST = round(gfst$Fst, 5)
    # Return mean and weighted FST.
  )
  
  flist[[length(flist)+1]] <- fst_dat
  
  progress(i) # % based progress bar.
  
}

# Join pair-specific FST outputs.
popdiff_df <- do.call("rbind", flist)

# Turn outputs from above into a symmetrical pairwise distance matrix.
flip_mat <- \(df) {
  
  df_mat <-  as.data.frame.matrix(xtabs(meanFST ~ ., df[,c(1:3)]))
  (mrow <- colnames(df_mat)[!colnames(df_mat) %in% rownames(df_mat)])
  (mcol <- rownames(df_mat)[!rownames(df_mat) %in% colnames(df_mat)])

  adjmat <- rbind(cbind(data.frame(mcol = rep(0, nrow(df_mat))), df_mat) %>% 
                       `colnames<-`(., c(mcol, colnames(df_mat))),
                  as.data.frame(matrix(data = 0, nrow = 1, ncol = ncol(df_mat)+1)) %>% 
                    `colnames<-`(., c(mcol, colnames(df_mat))) %>% 
                    `rownames<-`(., mrow))
                  
  adjmat[lower.tri(adjmat)] <- t(adjmat)[lower.tri(adjmat)]
  return(adjmat)
} 

# Pairwise distance matrix.
(pop_diff_mat <- flip_mat(popdiff_df))

# Write file and visualize.
write.csv(pop_diff_mat, "../data/hdgbs_maf005_fst_matrix.csv")
png(width = 1200, height = 1200, units = "px", "heatmap.png", res = 100)
gplots::heatmap.2(as.matrix(pop_diff_mat), trace = "none"); dev.off()




