#!/usr/bin/env Rscript

args <- commandArgs(TRUE)

# Confirm correct usage --------------------------------------------------------

if (length(args) != 2)
  stop("Usage: fst_pairwise.R <vcf> <site_info>")


# Install and load necessary packages ------------------------------------------

for (p in c("tidyverse", "SNPRelate", "gplots")) {
  if (!suppressMessages(require(p, character.only = T))) {
    message(paste("Installing:", p))
    install.packages(p, repos = "https://mirror.its.dal.ca/cran", dependencies = T)
    suppressMessages(require(p, character.only = T))}
  rm(p)
}

# Calculate Fst ---------------------------------------------------------------

# Close GDS envirnoment (if not starting fresh, this can cause issues with writing).
gdsfmt::showfile.gds(closeall = T)

# Name for output variables that are not easily overwritten.
name <- basename(gsub(pattern = "\\.vcf.*", replacement = "", args[1]))
print(name)

# Convert VCF to GDS object.
snpgdsVCF2GDS(vcf.fn = args[1],
              out.fn = `name`)

# Read in GDS object.
(genofile <- snpgdsOpen(`name`))
print("Genotype file read in correctly!")

# Read in individual sample information.
indvs <- read.csv(args[2], row.names = 1)

# Create a dataframe of all possible combinations of populations.
pop_pairs <- t(as.data.frame(combn(unique(indvs$site_full), 2)))

# List for populating.
flist <- list()

print("Calculating Fst")
# Calculates FST between each pair of populations.
for (i in 1:nrow(pop_pairs)) {
  
  pops <- pop_pairs[i,]
  
  # Samples of interest and populations (must be the same length).
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
  
}

print("Fst calculations complete")

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
pop_diff_mat <- flip_mat(df = popdiff_df)

# Write file and visualize.
write.csv(pop_diff_mat, paste0(name, ".csv"))
png(width = 1200, height = 1200, units = "px", paste0(name, "_heatmap.png"), res = 100)
gplots::heatmap.2(as.matrix(pop_diff_mat), trace = "none"); dev.off()
