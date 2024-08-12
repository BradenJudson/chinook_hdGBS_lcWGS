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
  mutate(Pop  = gsub("\\_.*","", basename(FileName)),
         MajAF = as.numeric(sub(".*:", "", MajAll)),
         MajNuc  = as.factor(substr(MajAll, start = 0, stop = 1)),
         MinAF = as.numeric(sub(".*:", "", MinAll)),
         MinNuc = as.factor(substr(MinAll, start = 0, stop = 1)),
         Gpos   = paste0(CHROM, "_", POS)) %>% as.data.frame() %>% 
  dplyr::select(c("Gpos", "MajAF")) %>% 
  pivot_wider(values_from = MajAF, names_from = Gpos)

temp  <- paste0("../data/pop_frq/6b_imputed_lcWGS_subset_134k/", 
                list.files(path = "../data/pop_frq/6b_imputed_lcWGS_subset_134k/"))
freqs <- as.data.frame(do.call(rbind, lapply(temp, rf))) %>% 
  mutate(Population = gsub("\\_.*","", basename(temp))) %>% 
  relocate(Population, 1)

write.table(freqs, "../data/pop_frq/6b_imputed_lcWGS_subset_134k.txt", 
            quote = F, row.names = F, col.names = T)


rmpops <- c("Big", "Harrison", "Raft", "Harrison")

freqs <- freqs %>% filter(!rownames(.) %in% rmpops)
  
# Site data --------------------------------------------------------------------

# Read in site information and format accordingly.
sites <- read.delim(file = "../data/ch2023_sequenced.txt") %>%
  mutate(Site = gsub("[[:space:]]", "", Site)) %>%
  relocate(Site, .after = last_col()) %>% 
  .[, c((ncol(.)-3):ncol(.))] %>% 
  filter(Site %in% rownames(freqs)) %>%
  arrange(Site)

sum(sites$Site == rownames(freqs)) == nrow(sites)

# Compute a linear distance (m) matrix accounting for the curvature of the earth.
# dm <- st_distance(st_as_sf(sites[,c(4,3,2)], coords = 2:3, crs = 4326))


rivdist <- read.csv("../data/Otsh_distances_mat.csv", row.names = 1) %>% 
  `colnames<-`(., gsub("\\.", "", colnames(.))) %>% 
  `rownames<-`(., gsub(" ", "", rownames(.))) 

rivdist <- data.matrix(rivdist)
rivdist <- rivdist[!rownames(rivdist) %in% c("Harrison", "Raft", "Nahatlatch"),
                   !colnames(rivdist) %in% c("Harrison", "Raft", "Nahatlatch")]

pc <- pcnm(dis = as.dist(rivdist))
pcnms <- as.data.frame(pc$vectors)

bioclim_sub <- as.matrix(cbind(subset(bioclim, select = bioc_to_use[,1]),
                               subset(pcnms, select = optpcnm[,1])))
                         



nmod  <- rda(freqs ~ 1, pcnms) # Null model.
fmod  <- rda(freqs ~ ., pcnms) # Full model.
optdb <- ordiR2step(nmod, scope = formula(fmod), direction = 'forward')
summary(optdb)$call # Prints formula with selected predictor variables.

png("../plots/pcnm_scores.png", width = 800, height = 1200, res = 100)
par(mfrow = c(3, 1))
hist(pcnms$PCNM1, main = "PCNM1", xlab = "PCNM Score")
hist(pcnms$PCNM3, main = "PCNM3", xlab = "PCNM Score")
hist(pcnms$PCNM4, main = "PCNM4", xlab = "PCNM Score")
dev.off()
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

sum(rownames(bioclim) %in% sites$Site) == nrow(sites)

# Heatmap/clustering of bioclimatic data and Latitude.
png("../plots/predictor_cors.png", width = 1200, height = 1000, res = 100)
bioclim_sub <- as.matrix(cbind(subset(bioclim, select = c( bio1, 
                                      bio5, bio8, bio16, bio15)),
                               subset(pcnms[,c(1,3:4)])))
gplots::heatmap.2(cor(bioclim_sub)^2, trace = 'none', margins = c(7,7)) ; dev.off()

sum(rownames(freqs) == rownames(bioclim)) == nrow(bioclim) 


################################################################################
# RDAs -------------------------------------------------------------------------
################################################################################

# Part I: hdGBS data -----------------------------------------------------------


# Define and run RDA.
(hdRDA <- vegan::rda(freqs ~ bio1 + bio5 + bio8 + bio16 + bio15 + Condition(PCNM1 + PCNM3 + PCNM4),
                     data = cbind(bioclim, pcnms), scale = T)) # Scale and center predictors.
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
  loadings <- scores(model, choices = c(1:4), display = 'species')
  png(paste0("../plots/", dataset, "_rda_loadings.png"), width = 1200, height = 1200, res = 100)
  par(mfrow = c(2, 2)); 
  hist(loadings[,1], main = "RDA1 Loadings", xlab = "RDA1")
  hist(loadings[,2], main = "RDA2 Loadings", xlab = "RDA2")
  hist(loadings[,3], main = "RDA3 Loadings", xlab = "RDA3")
  hist(loadings[,4], main = "RDA4 Loadings", xlab = "RDA4"); dev.off()
  print(paste0("Writing ", dataset, "_rda_loadings.png to ../plots/"), quote = F)
  
  # Identify and enumerate outlier loci along the first two RDA axes.
  cand1 <- outliers(loadings[,1], z = 3); print(paste0(length(cand1), " outliers on RDA1"), quote = F)
  cand2 <- outliers(loadings[,2], z = 3); print(paste0(length(cand2), " outliers on RDA2"), quote = F)
  cand3 <- outliers(loadings[,3], z = 3); print(paste0(length(cand3), " outliers on RDA3"), quote = F)
  cand4 <- outliers(loadings[,4], z = 3); print(paste0(length(cand4), " outliers on RDA4"), quote = F)
  
  # Isolate outlier loci into a pair of dataframes and adjust positional formatting.
  cand1DF <- cbind.data.frame(rep("RDA1", times=length(cand1)), names(cand1), unname(cand1))
  cand2DF <- cbind.data.frame(rep("RDA2", times=length(cand2)), names(cand2), unname(cand2))
  cand3DF <- cbind.data.frame(rep("RDA3", times=length(cand3)), names(cand3), unname(cand3))
  cand4DF <- cbind.data.frame(rep("RDA4", times=length(cand4)), names(cand4), unname(cand4))
  
  colnames(cand1DF) <- colnames(cand2DF) <- colnames(cand3DF) <- colnames(cand4DF) <- c("axis","snp","loading")
  
  all_cands <- rbind(cand1DF, cand2DF, cand3DF, cand4DF) %>% 
    separate(col = "snp", sep = 12,
             into = c("chromosome", "position"), remove = FALSE) %>% 
    mutate(chromosome = as.factor(gsub("_", "", chromosome)))
  
  # Only retain predictors (not conditioning variables).
  bio_vars_sub <- noquote(attr(rdamodel$terms, "term.labels"))[1:length(attr(rdamodel$terms, "term.labels"))-1]
  biodf <- bioclim[,c(bio_vars_sub)] # and the data itself
  
  # Make an empty matrix to populate in the following for loop.
  cand_mat <- matrix(nrow = nrow(all_cands), ncol = length(bio_vars_sub)) %>% 
    `colnames<-`(., c(bio_vars_sub))
  
  # Estimate the correlation between the allele frequency of each outlier SNP and each bioclim variable.
  # Important to make this an absolute value!
  for (i in 1:nrow(cand_mat)) {
    nam <- all_cands[i, 2]; snp <- freqs[,nam]
    cand_mat[i,] <- apply(biodf, 2, function(x) abs(cor(x, snp))) 
  } 
  
  # Add candidate SNP information to the correlation matrix computed above.
  cand_df_cor <- cbind.data.frame(all_cands, cand_mat) %>% rowwise()  %>% 
    mutate(correlation = max(c_across(starts_with("bio")), na.rm = TRUE), # Maximum correlation.
           correlation2 = sort(c_across(starts_with("bio")), TRUE)[[2]],  # Second highest correlation.
           predictor = names(.[,as.vector(bio_vars_sub)])[which.max(c_across(as.vector(bio_vars_sub)))],
           delta = abs(correlation - correlation2))  # Difference between max and second highest correlation.

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
  perc_var <- round(100*summary(rdamodel)$cont$importance[2, 1:4], 2)
  
  # Order site scores with respect to Latitude and join with RDA loadings.
  site_scores <- merge(sites, as.data.frame(scores(rdamodel, display = 'sites', scaling = 3, choices = c(1:4))),
                       by.x = "Site", by.y = 0) %>% arrange(Latitude) %>% mutate(Latitude = as.factor(Latitude))
  
  # Script for scaling = 3 plot.
  (rda_scale3_axes12 <- ggplot() +
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
  
  (rda_scale3_axes34 <- ggplot() +
      geom_point(data = as.data.frame(scores(rdamodel, display = 'species', 
                                             choices = c(3,4), scaling = 3)),
                 aes(x = RDA3, y = RDA4), colour = 'grey70', shape = 21, fill = 'grey90') +
      geom_point(data = site_scores, aes(x = RDA3, y = RDA4, fill = factor(Latitude)), shape = 21, size = 2) +
      scale_fill_manual(values = c(viridis_pal(option = "D")(length(unique(site_scores$Site)))),
                        labels = levels(unique(site_scores$Site))) +
      geom_segment(data = scores(rdamodel, display = 'bp', choices = c(3,4), scaling = 3),
                   aes(xend = RDA3*25, yend = RDA4*25, x = 0, y = 0), colour = '#0868ac',
                   lineend = "round", arrow = arrow(length = unit(0.1, "inches"))) +
      geom_text(data = as.data.frame(scores(rdamodel, display = 'bp', choices = c(3,4), 
                                            scaling = 3)) %>% rownames_to_column("bio"), colour = '#0868ac',
                aes(x = 26.5*RDA3, y = 26.5*RDA4, label = bio)) +
      theme_bw() + guides(fill = guide_legend(nrow = 5, byrow = TRUE)) +
      theme(legend.position = "none", legend.title = element_blank()) +
      labs(x = paste0("RDA3 (", perc_var[3], "%)"),
           y = paste0("RDA4 (", perc_var[4], "%)")))

  
  # Script for scaling = 1 plot.
  (rda_scale1_axes12 <- ggplot() +
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
  
  (rda_scale1_axes34 <- ggplot() +
      geom_point(data = RDAdf[is.na(RDAdf$correlation),], 
                 aes(x = RDA3, y = RDA4), colour = 'grey70',
                 shape = 21, fill = "gray95") +
      geom_point(data = RDAdf[!is.na(RDAdf$correlation),],
                 aes(x = RDA3, y = RDA4, fill = predictor),  
                 shape = 21, size = 2, alpha = 3/5) +
      scale_fill_manual(values = c(viridis_pal(option = "H")(length(unique(RDAdf$predictor)))),
                        labels = levels(unique(RDAdf$predictor))) +
      labs(x = paste0("RDA3 (", perc_var[3], "%)"),
           y = paste0("RDA4 (", perc_var[4], "%)")) +
      geom_segment(data = scores(rdamodel, display = 'bp', choices = c(3,4), scaling = 1),
                   aes(xend = RDA3, yend = RDA4, x = 0, y = 0), colour = '#0868ac', size = 1,
                   lineend = "round", arrow = arrow(length = unit(0.1, "inches"))) +
      geom_text(data = as.data.frame(scores(rdamodel, display = 'bp', choices = c(3,4), 
                                            scaling = 1)) %>% rownames_to_column("bio"), colour = '#0868ac',
                aes(x = 1.06*RDA3, y = 1.06*RDA4, label = bio)) + theme_bw() +
      theme(legend.title = element_blank(), legend.position = c(0.08, 0.8), legend.background = element_blank()))
  
  # Return base data frame and two ggplot objects as a list.
  return(list(RDAdf = RDAdf, rda_scale1_axes12 = rda_scale1_axes12, rda_scale1_axes34 = rda_scale1_axes34,
              rda_scale3_axes12 = rda_scale3_axes12, rda_scale3_axes34 = rda_scale3_axes34))
  
}

# Run above function on the hdGBS RDA defined earlier.
hdgbs_rda <- rda.full(model = hdRDA, data = "hdGBS")

# Plot both scaling = 1 and scaling = 3 plot types for the first 4 axes.
(hdgbs_plots <- cowplot::plot_grid(plotlist = list(
  hdgbs_rda$rda_scale3_axes12, hdgbs_rda$rda_scale1_axes12,
  hdgbs_rda$rda_scale3_axes34, hdgbs_rda$rda_scale1_axes34), 
  ncol = 2, nrow = 2, labels = c("a)", "b)", "c)", "d)")))
ggsave("../plots/hdgbs_rdas13.tiff", width = 16, height = 16)

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
  pivot_wider(values_from = MinAF, names_from = gPosition) %>%  
  arrange(Population) %>% column_to_rownames("Population")
  
rownames(lcafs) == rownames(bioclim)

# Part III: Imputed lcWGS data -------------------------------------------------

# Data are LD-pruned.
temp <- paste0("../data/pop_frq/3a_lcWGS_pruned_maf005_2M/", 
        list.files(path = "../data/pop_frq/3a_lcWGS_pruned_maf005_2M"))
implc <- as.data.frame(do.call(rbind, lapply(temp, rf))) %>% 
  `rownames<-`(., gsub("\\_.*","", basename(temp)))

write.csv(implc, "../data/pop_frq/imputed_lcwgs_afs.csv", row.names = T)
implc <- read.csv("../data/pop_frq/imputed_lcwgs_afs.csv") %>% 
  filter(rownames(.) %in% sites$site)

rownames(implc) == rownames(db)

(lcImpRDA <- vegan::rda(implc ~ bio2 + bio3 + bio8 + bio10 + bio13 + bio15 + Condition(dbMEM.1 + dbMEM.2 + dbMEM.8),
                        data = cbind(bioclim, db), scale = T))


# Part IV: Imputed and subset lcWGS --------------------------------------------

temp <- paste0("../data/pop_frq/6a_imputed_lcWGS_pruned_subset_105k/",
               list.files(path = "../data/pop_frq/6a_imputed_lcWGS_pruned_subset_105k"))
implcsub <- as.data.frame(do.call(rbind, lapply(temp, rf))) %>% 
  `rownames<-`(., gsub("\\_.*","", basename(temp)))

write_tsv(implcsub, "../data/pop_frq/6a_imputed_lcWGS_pruned_subset_105k.tsv")


sum(rownames(implcsub) == rownames(bioclim))

(lcImpRDA <- vegan::rda(implcsub ~ bio4 + bio2 + bio5 + bio8 + bio16 + bio15 + Condition(dbMEM.1 + dbMEM.2 + dbMEM.8),
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

a<-lcImpRDA_output$RDAdf %>% filter(!is.na(correlation)) %>% mutate(gpos = paste0(chromosome, position))
b<-hdgbs_rda$RDAdf %>% filter(!is.na(correlation)) %>% mutate(gpos = paste0(chromosome, position))

ggVennDiagram(x = list(a = a$gpos, b = b$gpos), category.names = c("lcImp", "hdgbs"))
