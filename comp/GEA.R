setwd("~/ots_landscape_genetics/comp")

library(tidyverse); library(vegan); #library(vcfR)
library(data.table); library(psych); library(broom)

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

temp  <- paste0("../data/pop_frq/", list.files(path = "../data/pop_frq"))
freqs <- as.data.frame(do.call(rbind, lapply(temp, rf))) %>% 
  `rownames<-`(., gsub("\\_.*","", str_sub(temp, start = 17)))


# Site data --------------------------------------------------------------------


sites <- read.delim(file = "../data/ch2023_sequenced.txt") %>%
  mutate(site = tools::toTitleCase(tolower(gsub("\\_.*", "", Population)))) %>%
  .[, c((ncol(.)-3):ncol(.))] %>% filter(site %in% rownames(bioclim)) %>% arrange(site)

sum(sites$site == rownames(bioclim)) == nrow(sites)

plot(sites[,2:3])

sites2 <- sites[,c(4,3,2)]
dm <- st_distance(st_as_sf(sites2, coords = 2:3, crs = 4326))
db <- create.dbMEM.model(coord = sites2[,c(2:3)], nsites = nrow(sites2)) %>% 
  as.data.frame() %>% `rownames<-`(., c(sites2$site))

memrd <- rda(freqs ~ 1, db)
memrf <- rda(freqs ~ ., db)
selmo <- ordiR2step(memrd, scope = formula(memrf), direction = 'both')

(hdRDA <- vegan::rda(freqs ~ bio2 +  bio8 + bio10 + bio13 + bio15 + Condition(dbMEM.1 + dbMEM.2 + dbMEM.8),
                     data = cbind(bioclim, db), scale = T)) # Scale and center predictors.


# Bioclimatic data -------------------------------------------------------------

# Read in and format bioclimatic data for each site.
bioclim <- read.csv("../data/ch_bioclim.csv") %>% 
  filter(is.na(ssp)) %>% .[,1:(ncol(.)-2)] %>% 
  mutate(Site = tools::toTitleCase(tolower(gsub("\\_.*", "", Site)))) %>% 
  filter(Site %in% rownames(freqs)) %>% arrange(Site) %>% column_to_rownames("Site")

# Heatmap/clustering of bioclimatic data and Latitude.
png("../plots/predictor_cors.png", width = 1200, height = 1000, res = 100)
bioclim_sub <- as.matrix(cbind(subset(bioclim, select = c(bio2, bio7, bio8,
                                      bio10, bio13, bio15)), Lat = sites[,2:3]))
gplots::heatmap.2(cor(bioclim_sub)^2, trace = 'none') ; dev.off()

sum(rownames(freqs) == rownames(bioclim)) == nrow(bioclim) 



# PCAs -------------------------------------------------------------------------

pcas <- read.csv("../data/hdgbs_sub_pcascores.csv") %>% 
  merge(., read.csv("../data/landgen_chinook_indvs.csv"),
        by.x = "X", by.y = "fish_ID") %>% 
  mutate(site_full = gsub("[[:space:]]Brood", "", site_full)) %>% 
  filter(site_full %in% rownames(bioclim_sub)) %>% 
  group_by(site_full) %>% arrange(site_full) %>% 
  summarise(PC1 = mean(PC1),
            PC2 = mean(PC2),
            PC3 = mean(PC3))

sum(sites$site == pcas$site_full) == nrow(sites)


################################################################################
# RDAs -------------------------------------------------------------------------
################################################################################

# Part I: hdGBS data -----------------------------------------------------------


# Define and run RDA.
(hdRDA <- vegan::rda(freqs ~ bio2 +  bio8 + bio10 + bio13 + bio15 + Condition(Latitude + Longitude),
                   data = cbind(bioclim, sites[,2:3]), scale = T)) # Scale and center predictors.
RsquareAdj(hdRDA) # Adjusted and non-adjusted model R2 values = 0.239 and 0.336, respectively.
(hd_terms <- anova.cca(hd, by = "terms")) # Reports significance of each term.
(hd_full  <- anova.cca(hd)) # Test and report significance of the full model.
(hd_axis  <- anova.cca(hd, by = "axis")) # Takes an incredible amount of time to run.
(eigvsumm <- round(summary(eigenvals(hdRDA, model = "constrained")), 3)) # % var explained by each axis.
screeplot(hdRDA) # Clearly indicates that RDA1 explains vast majority.
vif.cca(hdRDA) # Check model variance inflation factors.

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
  cand_mat <- matrix(nrow = nrow(cand), ncol = ncol(bio_vars_sub)) %>% 
    `colnames<-`(., c(colnames(bio_vars_sub)))
  
  # Estimate the correlation between the allele frequency of each outlier SNP and each bioclim variable.
  for (i in 1:nrow(cand_mat)) {
    nam <- all_cands[i, 2]; snp <- freqs[,nam]
    cand_mat[i,] <- apply(bio_vars_sub, 2, function(x) cor(x, snp)) } 
  
  # Add candidate SNP information to the correlation matrix computed above.
  cand_df_cor <- cbind.data.frame(all_cands, cand_mat) %>% rowwise() %>% 
    mutate(correlation = max(c_across(starts_with("bio")), na.rm = TRUE),
           predictor = names(.[6:11])[which.max(c_across(6:11))])
  
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
    mutate(predictor = case_when(
      !is.na(correlation) ~ predictor,
      is.na(correlation) ~ "non-outlier"))
  
  # Percent variation explained by each RDA axis.
  perc_var <- round(100*summary(rdamodel)$cont$importance[2, 1:2], 2)
  
  # Order site scores with respect to Latitude and join with RDA loadings.
  site_scores <- merge(sites, as.data.frame(scores(rdamodel, display = 'sites', scaling = 3)),
                       by.x = "site", by.y = 0) %>% arrange(Latitude) %>% mutate(Latitude = as.factor(Latitude))
  
  (rda_scale3 <- ggplot() +
      geom_point(data = as.data.frame(scores(rdamodel, display = 'species', 
                        choices = c(1,2), scaling = 3)),
                 aes(x = RDA1, y = RDA2), colour = 'grey70', shape = 21, fill = 'grey90') +
      geom_point(data = site_scores, aes(x = RDA1, y = RDA2, fill = factor(Latitude)), shape = 21, size = 2) +
      scale_fill_manual(values = c(viridis_pal(option = "D")(length(unique(site_scores$site)))),
                        labels = levels(unique(site_scores$site))) +
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
  
  
  return(list(RDAdf = RDAdf, rda_scale3 = rda_scale3))
  
}

hdgbs_rda <- rda.full(model = hdRDA, data = "hdGBS")
hdgbs_rda$rda_scale3; ggsave("../plots/hdgbs_rda_s3.tiff", width = 8, height = 8)


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
  



