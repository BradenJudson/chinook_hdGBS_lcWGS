setwd("~/ots_landscape_genetics/comp")

library(tidyverse); library(vegan); #library(vcfR)
library(data.table); library(psych); library(broom)

# From Brenna R. Forester at: https://popgen.nescent.org/2018-03-27_RDA_GEA.html
source("../scripts/outliers.R")

set.seed(240)


# TO DO -------------------------------------------------------------------

# 999 permutation sufficient or need to increase?


# HDBGS ------------------------------------------------------------------------


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


# Bioclimatic data -------------------------------------------------------------

# Read in and format bioclimatic data for each site.
bioclim <- read.csv("../data/ch_bioclim.csv") %>% 
  filter(is.na(ssp)) %>% .[,1:(ncol(.)-2)] %>% 
  mutate(Site = tools::toTitleCase(tolower(gsub("\\_.*", "", Site)))) %>% 
  filter(Site %in% rownames(freqs)) %>% arrange(Site) %>% column_to_rownames("Site")

# Heatmap/clustering of bioclimatic data and Latitude.
png("../plots/predictor_cors.png", width = 1200, height = 1000, res = 100)
bioclim_sub <- as.matrix(cbind(subset(bioclim, select = c(bio2, bio7, bio8,
                                      bio10, bio13, bio15)), Lat = sites$Latitude))
gplots::heatmap.2(cor(bioclim_sub)^2, trace = 'none') ; dev.off()

sum(rownames(freqs) == rownames(bioclim)) == nrow(bioclim) 

################################################################################
# RDAs -------------------------------------------------------------------------
################################################################################

# Part I: hdGBS data -----------------------------------------------------------


# Define and run RDA.
(hdRDA <- vegan::rda(freqs ~ bio2 + bio7 + bio8 + bio10 + bio13 + bio15, Z = Latitude, 
                   data = cbind(bioclim, sites[,2]), scale = T)) # Scale and center predictors.
RsquareAdj(hdRDA) # Adjusted and non-adjusted model R2 values = 0.239 and 0.336, respectively.
(hd_terms <- anova.cca(hd, by = "terms")) # Reports significance of each term.
(hd_full  <- anova.cca(hd)) # Test and report significance of the full model.
(hd_axis  <- anova.cca(hd, by = "axis")) # Takes an incredible amount of time to run.
(eigvsumm <- round(summary(eigenvals(hdRDA, model = "constrained")), 3)) # % var explained by each axis.
screeplot(hdRDA) # Clearly indicates that RDA1 explains vast majority.


### Design ggplot RDA function here ###

(hd_lds <- scores(hdRDA, choices = c(1:3), display = 'species'))
png("../plots/hd_rda_loadings.png", width = 1000, height = 600, res = 100)
par(mfrow = c(1, 2)); hist(hd_lds[,1], main = "RDA1 Loadings", xlab = "RDA1")
hist(hd_lds[,2], main = "RDA2 Loadings", xlab = "RDA2"); dev.off()

hd_cand1 <- outliers(hd_lds[,1], z = 3)
hd_cand2 <- outliers(hd_lds[,2], z = 3)

hd_cand1 <- cbind.data.frame(rep("Axis1", times=length(hd_cand1)), names(hd_cand1), unname(hd_cand1))
hd_cand2 <- cbind.data.frame(rep("Axis2", times=length(hd_cand2)), names(hd_cand2), unname(hd_cand2))
colnames(hd_cand1) <- colnames(hd_cand2) <- c("axis","snp","loading")
cand <- rbind(hd_cand1, hd_cand2) %>% 
  separate(col = "snp", sep = 12, into = c("chromosome", "position"),
           remove = FALSE) %>% 
  mutate(chromosome = as.factor(gsub("_", "", chromosome)))

# Remove Latitude as that is a conditioning variable in our RDA model.
bio_vars_sub <- as.data.frame(bioclim_sub) %>% select(starts_with("bio"))

pm <- matrix(nrow = nrow(cand), ncol = ncol(bio_vars_sub)) %>% 
  `colnames<-`(., c(colnames(bio_vars_sub)))

for (i in 1:nrow(cand)) {
  nam <- cand[i, 2]
  snp <- freqs[,nam]
  pm[i,] <- apply(bio_vars_sub, 2, function(x) cor(x, snp))
}

candDF <- cbind.data.frame(cand, pm)

for (i in 1:nrow(candDF)) {
  bar <- candDF[i,]
  candDF[i, 12] <- names(which.max(abs(bar %>% select(starts_with("bio")))))
  candDF[i, 13] <- max(abs(bar %>% select(starts_with("bio"))))
}; colnames(candDF)[12] <- 'predictor'; colnames(candDF)[13] <- 'correlation'

table(candDF$predictor)

dfwf <- candDF %>% pivot_wider(names_from = "axis", values_from = "loading") %>% 
  merge(hd_lds, by.x = "snp", by.y = 0)

j <- as.data.frame(hd_lds) %>% 
  merge(., dfwf[,c("snp", "predictor", "correlation")], 
        by.x = 0, by.y = "snp", all.x = TRUE) %>% 
  mutate(predictor = case_when(
    !is.na(correlation) ~ predictor,
     is.na(correlation) ~ "non-outlier"
  ))

ggplot(data = j, aes(x = RDA1, y = RDA2, fill = predictor)) + 
  geom_point(alpha= 1/5, shape = 21)


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
  



