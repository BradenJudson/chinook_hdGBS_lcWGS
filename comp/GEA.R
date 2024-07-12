setwd("~/ots_landscape_genetics/comp")

library(tidyverse); library(vegan); #library(vcfR)
library(data.table); library(psych); library(broom)

set.seed(240)

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
freqs <- do.call(rbind, lapply(temp, rf)) %>% 
  `rownames<-`(., gsub("\\_.*","", str_sub(temp, start = 17)))


# Bioclimatic data -------------------------------------------------------------


bioclim <- read.csv("../data/ch_bioclim.csv") %>% 
  filter(is.na(ssp)) %>% .[,1:(ncol(.)-2)] %>% 
  mutate(Site = tools::toTitleCase(tolower(gsub("\\_.*", "", Site)))) %>% 
  filter(Site %in% rownames(freqs)) %>% arrange(Site) %>% column_to_rownames("Site")

as.data.frame(cor(bioclim) %>% replace(. < 0.8, NA)); pairs.panels(bioclim, scale = T)


bioclim_sub <- subset(bioclim, select = -c(bio11, bio4, bio6, bio7, bio9, bio13, 
                                           bio12, bio14, bio19, bio17, bio10,
                                           bio16))

pairs.panels(bioclim_sub)

sum(rownames(freqs) == rownames(bioclim)) == nrow(bioclim) 

# Below is trying and failing to randomly select uncorrelated variables from a 
#correlation matrix.
# Idea - collapse everything long form, get all unique combinations of variables, 
# group by and randomly subset to 1?
# 
# cb <- abs(cor(bioclim))
# cb[lower.tri(cb, diag = T)] <- NA
# 
# jp <- list()
# for (i in colnames(cb)) {
#   cv <- c(cb[,i] > 0.8)
#   da <- data.frame(v1 = names(cb[cv, i][1:sum(cv, na.rm = T)]),
#                    v2 = i) %>% 
#     group_by(v2) %>% 
#     dplyr::sample_n(1)
#   
#   jp[[length(jp) + 1]] <- da
#   
#     
# }
# 
# k <- do.call(rbind, jp) %>% filter(!is.na(v1))
# 

# Site data --------------------------------------------------------------------

sites <- read.delim(file = "../data/ch2023_sequenced.txt") %>%
  mutate(site = tools::toTitleCase(tolower(gsub("\\_.*", "", Population)))) %>%
  .[, c((ncol(.)-3):ncol(.))] %>% filter(site %in% rownames(bioclim)) %>% arrange(site)

sum(sites$site == rownames(bioclim)) == nrow(sites)


# Population structure ---------------------------------------------------------

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


# RDAs -------------------------------------------------------------------------

# Define possible predictor variables.
preds <- cbind(bioclim_sub, pcas) %>% 
  select(-c("site_full"))

# Run models as defined and assign their respective ANOVA to a similarly named object.
(m1 <- vegan::rda(freqs ~ ., data = preds[,c(1:5, 10:11)], scale = T))
(m1_aov <- anova.cca(m1, by = "terms", parallel=getOption("mc.cores"))); RsquareAdj(m1)

(m2 <- vegan::rda(freqs ~ ., Z = PC1, data = preds[,c(2:8)], scale = T))
(m2_aov <- anova.cca(m2, by = "terms", parallel=getOption("mc.cores"))); RsquareAdj(m2)

(m3 <- vegan::rda(freqs ~ bio2 + bio5 + bio8 + bio15 + bio18, Z = Latitude, data = preds, scale = T))
(m3_aov <- anova.cca(m3, by = "terms", parallel=getOption("mc.cores"))); RsquareAdj(m3)

(m4 <- vegan::rda(freqs ~ bio2 + bio5 + bio8 + bio15 + bio18 + Condition(PC1 + PC2), data = preds, scale = T))
(m4_aov <- anova.cca(m4, by = "terms", parallel=getOption("mc.cores"))); RsquareAdj(m4)

(m5 <- vegan::rda(freqs ~ bio2 + bio5 + bio8 + bio15 + bio18 + Condition(Latitude + Longitude), data = preds, scale = T))
(m5_aov <- anova.cca(m5, by = "terms", parallel=getOption("mc.cores"))); RsquareAdj(m5)

(m6 <- vegan::rda(freqs ~ bio2 + bio5 + bio8 + bio15 + bio18 + Condition(PC1 + PC2 + PC3), data = preds, scale = T))
(m6_aov <- anova.cca(m6, by = "terms", parallel=getOption("mc.cores"))); RsquareAdj(m6)

(m7 <- vegan::rda(freqs ~ ., Z = PC1, data = preds[,c(2:8,10)], scale = T))
(m7_aov <- anova.cca(m7, by = "terms", parallel=getOption("mc.cores"))); RsquareAdj(m7)

(m8 <- vegan::rda(freqs ~ bio2 + bio3 + bio5 + bio8 + bio15 + bio16 + bio18 + Condition(PC1 + PC2 + PC3), data = preds, scale = T))
(m8_aov <- anova.cca(m8, by = "terms")); RsquareAdj(m8)


# m0 <- rda(freqs ~ Condition(PC1 + PC2 + PC3), data = preds)
# m1 <- rda(freqs ~ Condition(PC1 + PC2 + PC3) + bio1 + bio2 + bio3 + bio5 + bio8 + bio15 + bio18, data = preds)

(frda <- rda(freqs ~ bio1 + bio2 + bio3 + bio4 + bio5 + bio6 + bio7 + bio8 + bio9 + bio10 + bio11 + bio12 + bio13 + bio14 + bio15 + bio16 + bio17 + bio18 + bio19 + Condition(PC1 + PC2 + PC3), 
            data = cbind(bioclim, pcas[,c(2:4)]), scale = T))

(fwdsel <- ordiR2step(rda(freqs ~ 1 + Condition(PC1 + PC2 + PC3), 
                          data = cbind(bioclim, pcas[,c(2:4)])),
                      scope = formula(frda),
                      direction = "forward",
                      R2scope = TRUE,
                      pstep = 1000,
                      trace = TRUE))

(m9 <- vegan::rda(freqs ~ bio1 + bio15 + bio7 + bio9, data = bioclim, scale = T))
(m9_aov <- anova.cca(m9)); RsquareAdj(m9)

(m10 <- vegan::rda(freqs ~ bio1 + bio15 + bio7 + bio9 + Condition(PC1 + PC2 + PC3), 
                   data = cbind(bioclim, pcas[,c(2:4)]), scale = T))
(m10_aov <- anova.cca(m10, by = "terms")); RsquareAdj(m10)


(m11 <- vegan::rda(freqs ~ bio19 + bio14 + bio5, data = cbind(bioclim, pcas[,c(2:4)]), scale = T))
(m11_aov <- anova.cca(m11, by = "terms")); RsquareAdj(m11)

(m12 <- vegan::rda(freqs ~ bio1 + bio15 + bio7 + bio9 + Condition(PC1 + PC2), 
                   data = cbind(bioclim, pcas[,c(2:4)]), scale = T))
(m12_aov <- anova.cca(m12)); RsquareAdj(m12)


rda_tidy <- \(x) { # Function for summarizing RDA models.

  # Retrieve RDA and corresponding ANOVA objects.
  rda_name  <- rlang::as_label(eval(parse(text=enquo(x)))[[2]])
  rda_model <- get(x = rda_name, envir = .GlobalEnv)
  aov       <- get(x = paste0(rda_name, "_aov"), envir = .GlobalEnv)
  aov_stats <- broom::tidy(aov)[1, c("df", "statistic", "p.value")]
  formula   <- gsub(", data = preds, scale = T","", 
                    str_sub(string = gsub(")","", 
                    c(rda_model$call)), start = 15))
  AdjR2     <- as.numeric(RsquareAdj(rda_model)[2])

  # Combine summary stats and equation into a dataframe.
  df <- data.frame(m = rda_name,
  eq = formula, AdjR2 = AdjR2
  ) %>% cbind(., aov_stats)
  
  return(df)

}

# Return RDA summaries in a single dataframe.
rda_summaries <- noquote(paste0("m", seq(1,2))) %>% 
  lapply(., rda_tidy) %>% do.call(rbind, .)



# ------------------------------------------------------------------------------

# Retrieve genotype matrix from vcf [THIS HAS MISSING DATA - NEED TO FIX]
hdgbs <- read.vcfR("../data/snps_maf001_singletons.vcf")
gts   <- t(extract.gt(hdgbs, element = "GT", as.numeric = TRUE)) %>% 
  as.data.frame() %>% `rownames<-`(., gsub("\\_.*", "", rownames(.)))



