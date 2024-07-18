setwd("~/ots_landscape_genetics/comp")

library(tidyverse); library(vegan); library(sf); library(adespatial)
library(data.table); library(psych); library(viridis)

# From Brenna R. Forester at: https://popgen.nescent.org/2018-03-27_RDA_GEA.html
source("../scripts/outliers.R")

set.seed(240) 


# TO DO -------------------------------------------------------------------

# 999 permutation sufficient or need to increase?


# HDBGS ------------------------------------------------------------------------

# Function for reading and formatting the output of vcftools --freq.
rf <- function(dir) paste0(dir, list.files(pattern = ".*.frq", 
                           path = dir)) %>% set_names(.) %>% 
  map_df(read_table, .id = "FileName", 
         col_names = c("CHROM", "POS", "N_ALLELES", 
                       "N_CHR", "MajAll", "MinAll")) %>% 
  filter(MajAll != "{ALLELE:FREQ}") %>% 
  mutate(Pop  = gsub("\\_.*","", str_sub(FileName, start = 17)),
         MajAF = as.numeric(sub(".*:", "", MajAll)),
         MajNuc  = as.factor(substr(MajAll, start = 0, stop = 1)),
         MinAF = as.numeric(sub(".*:", "", MinAll)),
         MinNuc = as.factor(substr(MinAll, start = 0, stop = 1)),
         Gpos   = paste0(CHROM, "_", POS)) %>% 
  select(c("Gpos", "MajAF")) %>% 
  pivot_wider(values_from = MajAF, names_from = Gpos)

temp  <- paste0("../data/pop_frq/hdgbs_not_pruned/", 
                list.files(path = "../data/pop_frq/hdgbs_not_pruned"))
freqs <- as.data.frame(do.call(rbind, lapply(temp, rf))) %>% 
  `rownames<-`(., gsub("\\_.*","", str_sub(temp, start = 34)))


# Site data --------------------------------------------------------------------

# Read in site information and format accordingly.
sites <- read.delim(file = "../data/ch2023_sequenced.txt") %>%
  mutate(Site = gsub("[[:space:]]", "", Site)) %>%
  relocate(Site, .after = last_col()) %>% 
  .[, c((ncol(.)-3):ncol(.))] %>% filter(Site %in% rownames(freqs)) %>% arrange(Site)

sum(sites$Site == rownames(freqs)) == nrow(sites)


# Compute a linear distance (m) matrix accounting for the curvature of the earth.
dm <- st_distance(st_as_sf(sites[,c(4,3,2)], coords = 2:3, crs = 4326))

# Convert distance matrix into distance-based Moran's eigenvector maps.
db <- create.dbMEM.model(coord = sites[,c(2:3)], nsites = nrow(sites)) %>% 
  as.data.frame() %>% `rownames<-`(., c(sites$site)) # Formatting.

nmod  <- rda(freqs ~ 1, db) # Null model.
fmod  <- rda(freqs ~ ., db) # Full model.
optdb <- ordiR2step(nmod, scope = formula(fmod), direction = 'both')
summary(optdb)$call # Prints formula with selected predictor variables.
# Choose dbMEMs that optimally describe spatial variation in response data.

# Bioclimatic data -------------------------------------------------------------

# Read in and format bioclimatic data for each site.
bioclim <- read.csv("../data/ch_bioclim.csv") %>% 
  filter(is.na(ssp)) %>% .[,1:(ncol(.)-2)] %>% 
  mutate(Site = gsub("Qualicum", "BigQualicum", gsub(" ", "", 
                     tools::toTitleCase(tolower(gsub("\\_", " ", 
                gsub("_FALL", "", gsub("_HATCHERY_FALL" ,"", 
                gsub("\\-.*", "", gsub("\\_CREEK*", "", 
                gsub("\\_RIVER*", "", Site))))))))))) %>% 
  filter(Site %in% rownames(freqs)) %>% arrange(Site) %>% 
  column_to_rownames("Site")

sum(bioclim$Site %in% sites$Site) == nrow(sites)

# Heatmap/clustering of bioclimatic data and Latitude.
png("../plots/predictor_cors.png", width = 1200, height = 1000, res = 100)
bioclim_sub <- as.matrix(cbind(subset(bioclim, select = c(bio2, bio1, 
                                      bio10, bio13, bio15)),
                               subset(db, select = c(dbMEM.1, dbMEM.2, dbMEM.8))))
gplots::heatmap.2(cor(bioclim_sub)^2, trace = 'none', margins = c(7,7)) ; dev.off()

sum(rownames(freqs) == rownames(bioclim)) == nrow(bioclim) 


################################################################################
# RDAs -------------------------------------------------------------------------
################################################################################

# Part I: hdGBS data -----------------------------------------------------------


# Define and run RDA.
(hdRDA <- vegan::rda(freqs ~ bio4 + bio2 + bio5 + bio8 + bio16 + bio15 + Condition(dbMEM.1 + dbMEM.2 + dbMEM.8),
                     data = cbind(bioclim, db), scale = T)) # Scale and center predictors.
RsquareAdj(hdRDA) # Adjusted and non-adjusted model R2 values = 0.158 and 0.067, respectively.
(hd_terms <- anova.cca(hdRDA, by = "terms")) # Reports significance of each term.
(hd_full  <- anova.cca(hdRDA)) # Test and report significance of the full model.
(hd_axis  <- anova.cca(hdRDA, by = "axis"))  # Takes an incredible amount of time to run. RDA1 and 2 significant.
(eigvsumm <- round(summary(eigenvals(hdRDA, model = "constrained")), 3)) # % var explained by each axis.
screeplot(hdRDA) # Variation explained by each RDA axis.
round(vif.cca(hdRDA), 2) # Check model variance inflation factors.

# Define a function for extracting RDA elements and plotting the results.
# Will use this same function with the other data sources too.
rda.full <- \(model, dataset) {
  
  # Pull in the RDA model of interest.
  rdamodel <- get(x = as_label(eval(parse(text=enquo(model)))[[2]]),
                  envir = .GlobalEnv)
  
  # SNP loadings across first two RDA axes and the distributions of those loadings. 
  loadings <- scores(model, choices = c(1:2), display = 'species')
  png(paste0("../plots/", dataset, "_rda_loadings.png"), width = 1400, height = 600, res = 100)
  par(mfrow = c(1, 2)); hist(loadings[,1], main = "RDA1 Loadings", xlab = "RDA1")
  hist(loadings[,2], main = "RDA2 Loadings", xlab = "RDA2"); dev.off()
  print(paste0("Writing ", dataset, "_rda_loadings.png to ../plots/"), quote = F)
  
  # Identify and enumerate outlier loci along the first two RDA axes.
  cand1 <- outliers(loadings[,1], z = 3); print(paste0(length(cand1), " outliers on RDA1"), quote = F)
  cand2 <- outliers(loadings[,2], z = 3); print(paste0(length(cand2), " outliers on RDA2"), quote = F)
  
  # Isolate outlier loci into a pair of dataframes and adjust positional formatting.
  cand1DF <- cbind.data.frame(rep("RDA1", times=length(cand1)), names(cand1), unname(cand1))
  cand2DF <- cbind.data.frame(rep("RDA2", times=length(cand2)), names(cand2), unname(cand2))
  colnames(cand1DF) <- colnames(cand2DF) <- c("axis","snp","loading")
  all_cands <- rbind(cand1DF, cand2DF) %>% 
    separate(col = "snp", sep = 12,
             into = c("chromosome", "position"), remove = FALSE) %>% 
    mutate(chromosome = as.factor(gsub("_", "", chromosome)))
  
  # Only retain predictors (not conditioning variables).
  bio_vars_sub <- as.data.frame(bioclim_sub) %>% select(starts_with("bio"))
  
  # Make an empty matrix to populate in the following for loop.
  cand_mat <- matrix(nrow = nrow(all_cands), ncol = ncol(bio_vars_sub)) %>% 
    `colnames<-`(., c(colnames(bio_vars_sub)))
  
  # Estimate the correlation between the allele frequency of each outlier SNP and each bioclim variable.
  for (i in 1:nrow(cand_mat)) {
    nam <- all_cands[i, 2]; snp <- freqs[,nam]
    cand_mat[i,] <- apply(bio_vars_sub, 2, function(x) cor(x, snp)) } 
  
  # Add candidate SNP information to the correlation matrix computed above.
  cand_df_cor <- cbind.data.frame(all_cands, cand_mat) %>% rowwise() %>% 
    mutate(correlation = max(c_across(starts_with("bio")), na.rm = TRUE),
           predictor = names(.[6:(5+ncol(bio_vars_sub))])[which.max(c_across(6:(5+ncol(bio_vars_sub))))])  

  # Count how many outlier SNPs associate with each predictor variable.
  print(table(cand_df_cor$predictor))
  
  # Pivot correlation matrix to long form and merge SNP loadings along each RDA axis.
  dfwf <- cand_df_cor %>% 
    pivot_wider(names_from = "axis", values_from = "loading") %>% 
    merge(loadings, by.x = "snp", by.y = 0)
  
  # Join outlier SNPs with non-outlier SNPs for plotting purposes.
  RDAdf <- as.data.frame(loadings) %>% 
    merge(., dfwf[,c("snp", "predictor", "correlation")], 
          by.x = 0, by.y = "snp", all.x = TRUE) %>% 
    mutate(predictor = as.factor(case_when(
      !is.na(correlation) ~ predictor,
      is.na(correlation) ~ "non-outlier"))) %>% 
    separate(col = Row.names, c(12), into = c("chromosome", "position")) %>% 
    mutate(chromosome = as.factor(gsub("_", "", chromosome)),
           position = as.numeric(position))
  
  # Percent variation explained by each RDA axis.
  perc_var <- round(100*summary(rdamodel)$cont$importance[2, 1:2], 2)
  
  # Order site scores with respect to Latitude and join with RDA loadings.
  site_scores <- merge(sites, as.data.frame(scores(rdamodel, display = 'sites', scaling = 3)),
                       by.x = "Site", by.y = 0) %>% arrange(Latitude) %>% mutate(Latitude = as.factor(Latitude))
  
  # Script for scaling = 3 plot.
  (rda_scale3 <- ggplot() +
      geom_point(data = as.data.frame(scores(rdamodel, display = 'species', 
                        choices = c(1,2), scaling = 3)),
                 aes(x = RDA1, y = RDA2), colour = 'grey70', shape = 21, fill = 'grey90') +
      geom_point(data = site_scores, aes(x = RDA1, y = RDA2, fill = factor(Latitude)), shape = 21, size = 2) +
      scale_fill_manual(values = c(viridis_pal(option = "D")(length(unique(site_scores$Site)))),
                        labels = levels(unique(site_scores$Site))) +
      geom_segment(data = scores(rdamodel, display = 'bp', choices = c(1,2), scaling = 3),
                 aes(xend = RDA1*25, yend = RDA2*25, x = 0, y = 0), colour = '#0868ac',
                 lineend = "round", arrow = arrow(length = unit(0.1, "inches"))) +
      geom_text(data = as.data.frame(scores(rdamodel, display = 'bp', choices = c(1,2), 
                       scaling = 3)) %>% rownames_to_column("bio"), colour = '#0868ac',
                aes(x = 26.5*RDA1, y = 26.5*RDA2, label = bio)) +
      theme_bw() + guides(fill = guide_legend(nrow = 5, byrow = TRUE)) +
      theme(legend.position = "none", legend.title = element_blank()) +
      labs(x = paste0("RDA1 (", perc_var[1], "%)"),
           y = paste0("RDA2 (", perc_var[2], "%)")))
  ggsave(paste0("../plots/", dataset, "_rda_s3.tiff"), width = 8, height = 8)
  print(paste0("Writing ", dataset, "_rda_s3.tiff to ../plots/"), quote = F)
  
  # Script for scaling = 1 plot.
  (rda_scale1 <- ggplot() +
      geom_point(data = RDAdf[is.na(RDAdf$correlation),], 
                 aes(x = RDA1, y = RDA2), colour = 'grey70',
                 shape = 21, fill = "gray95") +
      geom_point(data = RDAdf[!is.na(RDAdf$correlation),],
                 aes(x = RDA1, y = RDA2, fill = predictor),  
                 shape = 21, size = 2, alpha = 3/5) +
      scale_fill_manual(values = c(viridis_pal(option = "H")(length(unique(RDAdf$predictor)))),
                        labels = levels(unique(RDAdf$predictor))) +
      labs(x = paste0("RDA1 (", perc_var[1], "%)"),
           y = paste0("RDA2 (", perc_var[2], "%)")) +
      geom_segment(data = scores(rdamodel, display = 'bp', choices = c(1,2), scaling = 1),
                   aes(xend = RDA1, yend = RDA2, x = 0, y = 0), colour = '#0868ac', size = 1,
                   lineend = "round", arrow = arrow(length = unit(0.1, "inches"))) +
      geom_text(data = as.data.frame(scores(rdamodel, display = 'bp', choices = c(1,2), 
                                            scaling = 1)) %>% rownames_to_column("bio"), colour = '#0868ac',
                aes(x = 1.06*RDA1, y = 1.06*RDA2, label = bio)) + theme_bw() +
      theme(legend.title = element_blank(), legend.position = c(0.08, 0.8), legend.background = element_blank())) 
  
  # Return base data frame and two ggplot objects as a list.
  return(list(RDAdf = RDAdf, rda_scale3 = rda_scale3, rda_scale1 = rda_scale1))
  
}

# Run above function on the hdGBS RDA defined earlier.
hdgbs_rda <- rda.full(model = hdRDA, data = "hdGBS")

# Plot both scaling = 1 and scaling = 3 plot types. 
(hdgbs_plots <- cowplot::plot_grid(plotlist = list(
  hdgbs_rda$rda_scale3, hdgbs_rda$rda_scale1), 
  ncol = 2, labels = c("a)", "b)")))
ggsave("../plots/hdgbs_rdas13.tiff", width = 16, height = 8)

# Quick function to plot where the outliers are throughout the genome.
outlier_pos <- \(full_rda) {
  
  # Base plot of chromosome lengths, and a dataframe for the outliers only.
  genbar <- full_rda$RDAdf %>% group_by(chromosome) %>% summarise(chLen = max(position))
  hd_out <- full_rda$RDAdf %>% filter(!is.na(correlation))
  
  (plot <- ggplot(data = NULL) +
      geom_bar(data = genbar, 
               aes(x = chromosome, y = chLen/1e6),
               stat = 'identity', fill = 'gray90',
               colour = 'gray90', width = 1/5) +
      geom_point(data = hd_out, 
                 aes(x = chromosome, y = position/1e6)) +
      labs(y = "Base pairs (millions)", x = NULL) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, vjust = 1/2)))
}

(hdgbs_outliers_positions <- outlier_pos(hdgbs_rda))
ggsave("../plots/hdgbs_outlier_positions.tiff", width = 10, height = 6)



# Part II: lcWGS data -----------------------------------------------------

# Read in population-specific allele frequencies derived from genotype likelihoods.
# Requires some awkward formatting and transposition to prepare for RDAs.
lcafs <- paste0("../lcWGS/pop_afs/", list.files(pattern = "*.txt", 
                                path = "../lcWGS/pop_afs/")) %>% set_names(.) %>% 
  map_df(read_table, .id = "Population") %>% rename("MinAF" = knownEM) %>% 
  mutate(Population = tools::toTitleCase(gsub("\\.txt", "", str_sub(start = 18, Population))),
         gPosition  = paste0(chromo, "_", position)) %>% 
  select(c("Population", "MinAF", "gPosition")) %>% 
  pivot_wider(values_from = MinAF, names_from = gPosition) %>%     #### USE ARRANGE HERE TO GET SAME AS HDGBS ####
  column_to_rownames("Population")
  

# Part III: Imputed lcWGS data -------------------------------------------------

# Data are LD-pruned.
temp <- paste0("../data/pop_frq/imputed_lcwgs/", 
        list.files(path = "../data/pop_frq/imputed_lcwgs"))
implc <- as.data.frame(do.call(rbind, lapply(temp, rf))) %>% 
  `rownames<-`(., gsub("\\_.*","", str_sub(temp, start = 31)))

# write.csv(implc, "../data/pop_frq/imputed_lcwgs_afs.csv", row.names = T)
implc <- read.csv("../data/pop_frq/imputed_lcwgs_afs.csv") %>% 
  filter(rownames(.) %in% sites$site)

rownames(implc) == rownames(db)

(lcImpRDA <- vegan::rda(implc ~ bio2 + bio3 + bio8 + bio10 + bio13 + bio15 + Condition(dbMEM.1 + dbMEM.2 + dbMEM.8),
                        data = cbind(bioclim, db), scale = T))


# Part IV: Imputed and subset lcWGS --------------------------------------------

temp <- paste0("../data/pop_frq/imputed_hdgbs_subset/", 
               list.files(path = "../data/pop_frq/imputed_hdgbs_subset"))
implcsub <- as.data.frame(do.call(rbind, lapply(temp, rf))) %>% 
  `rownames<-`(., gsub("San", "SanJuan", gsub("\\_.*","", str_sub(temp, start = 38)))) %>% 
  filter(rownames(.) %in% rownames(bioclim)) %>% arrange(rownames(.))




sum(rownames(implcsub) == rownames(bioclim))

(lcImpRDA <- vegan::rda(implcsub ~ bio2 + bio8 + bio10 + bio13 + bio15 + Condition(dbMEM.1 + dbMEM.2 + dbMEM.8),
                        data = cbind(bioclim, db), scale = T))
vif.cca(lcImpRDA)
RsquareAdj(lcImpRDA)
lcImpRDA_output <- rda.full(model = lcImpRDA, dataset = "imputed_lcWGS")

(lcimpsub_plots <- cowplot::plot_grid(plotlist = list(
  lcImpRDA_output$rda_scale3, lcImpRDA_output$rda_scale1), 
  ncol = 2, labels = c("a)", "b)")))
ggsave("../plots/lcImpRDA_RDA13.tiff", width = 16, height = 8)

################################################################################
# Comparisons ------------------------------------------------------------------
################################################################################

# Outlier comparisons ----------------------------------------------------------

