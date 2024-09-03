setwd("~/ots_landscape_genetics/hdGBS/imputation_comparison")

library(vcfR); library(tidyverse)

# Part I: Replace missing with most common genotype ----------------------------

# Read in original VCF with missing data.
original <- read.vcfR("vcfs/hdgbs_subset_134kSNPs_n362_original.vcf")
o.mat <- extract.gt(original, return.alleles = TRUE) # Extract genotypes.

o.mat[o.mat == "."] <- NA
sum(is.na(o.mat)) # Number of missing sites.
round(sum(is.na(o.mat))/(nrow(o.mat)*ncol(o.mat))*100, 2) # Proportion of missing sites.

missing <- which(is.na(o.mat), T) # Find positions of all missing SNPs.

# Replace missing genotypes with most common genotype at that site.
mcg <- apply(o.mat, 2, function(x) replace(x, is.na(x), names(which.max(table(x))))) 

# Isolate imputed genotypes using this approach.
mcsub <- mcg[missing]

# Convert to allele frequency matrix.
snp_info <- original@fix

mcg_aflist <- list()

for (i in 1:ncol(mcg)) {
  
  df <- data.frame(genotypes = mcg[,i]) %>% 
    rownames_to_column("ID") %>% 
    merge(snp_info, by = "ID") %>% 
    mutate(A1 = str_sub(genotypes,1,1),
           A2 = str_sub(genotypes,3,3)) %>% 
    mutate(GPOS = paste0(CHROM, "_", POS),
           fish = colnames(mcg)[i])
  
  mcg_aflist[[length(mcg_aflist) + 1]] <- df
  
}

afs <- as.data.frame(do.call(rbind, mcg_aflist)) 

sum(afs$GPOS %in% dat$GPOS) == nrow(afs)

afsm <- merge(afs, dat[,c("MajNuc", "MinNuc", "GPOS")], by = "GPOS")

hd_mcgeno  <- cbind(afs %>% arrange(GPOS),
              dat[,c("MajNuc", "MinNuc", "GPOS")] %>% 
              arrange(GPOS)) %>% 
  select(c("ID", "genotypes", "CHR","POS", "A1", 
           "A2", "MajNuc", "MinNuc", "fish")) %>% 
  mutate(freq = case_when(A1 == MajNuc & A2 != MajNuc | A1 != MajNuc & A2 == MajNuc ~ 0.5,
                          A1 != MajNuc & A2 != MajNuc ~ 0.0,
                          A1 == MajNuc & A2 == MajNuc ~ 1.0)) %>% 
  select(c("GPOS", "fish", "freq")) %>% 
  pivot_wider(names_from = fish, 
              values_from = freq) %>% 
  group_by(GPOS) %>%  
  summarise(AF = mean(c_across(2:ncol(.)-1), na.rm = TRUE))


# Part II: Use BEAGLE to imputed missing genotypes -----------------------------

# Read in VCF where missing data were imputed using BEAGLE v5.
imputed <- read.vcfR("vcfs/hdgbs_subset_134kSNPs_n362_imputed.vcf")
i.mat <- extract.gt(imputed) # Extract genotypes.

# Isolate imputed genotypes only.
isub <- i.mat[missing]; length(isub)

# Confirm ordering and dimensions are consistent.
all.equal(nrow(mcg), nrow(o.mat), nrow(i.mat))
all.equal(length(missing)/2, length(mcsub), length(isub))

# Combine everything into a single dataframe.
# Requires some reformatting.
full <- as.data.frame(missing) %>% 
  mutate(beagle = sub("\\|", "/", isub),
         mostco = mcsub)

# Variation in imputation methods.
(match  <- sum(full$beagle == full$mostco))  
(nmatch <- sum(full$beagle != full$mostco))  
(total  <- nrow(full))                       
(prop   <- round(match/total*100, 2))        

# Comparisons ------------------------------------------------------------------

rf <- function(dir) paste0(dir, 
                           list.files(pattern = ".*.txt.frq", 
                                      path = dir)) %>% set_names(.) %>% 
  map_df(read_table, .id = "FileName", 
         col_names = c("CHROM", 
                       "POS", "N_ALLELES", 
                       "N_CHR", "MajAll", 
                       "MinAll")) %>% 
  filter(MajAll != "{ALLELE:FREQ}") %>% 
  mutate(sample = gsub("_freq.txt.frq", "", gsub(dir, "", FileName)),
         MajAF = as.numeric(sub(".*:", "", MajAll)),
         MajNuc  = as.factor(substr(MajAll, start = 0, stop = 1)),
         MinAF = as.numeric(sub(".*:", "", MinAll)),
         MinNuc = as.factor(substr(MinAll, start = 0, stop = 1))) %>% 
  select(!c("MajAll", "MinAll", "N_CHR", "N_ALLELES", "FileName"))

hdBeagle <- rf(dir = "./indv_afs_134k/hdgbs_beagle/")
lcBeagle <- rf(dir = "./indv_afs_134k/lcwgs_beagle/") %>% 
  rename(MajAFlc = "MajAF", MajNuclc = "MajNuc", 
         MinAFlc = "MinAF", MinNuclc = "MinNuc")

dat <- merge(hdBeagle, lcBeagle, by = c("CHROM", "POS", "sample"))
write.csv(dat, "hd_lc_comb_pre_fix.csv", row.names = F)

# Rearrange columns such that major and minor alleles are consistent.
dat[dat$MajNuc==dat$MinNuclc & dat$MinNuc==dat$MajNuclc,  ] <- dat[dat$MajNuc==dat$MinNuclc & dat$MinNuc==dat$MajNuclc,   c(1:3,6:7,4:5,8:11) ]
dat[(dat$MajNuc==dat$MinNuclc & dat$MinNuc!=dat$MajNuclc),] <- dat[(dat$MajNuc==dat$MinNuclc & dat$MinNuc!=dat$MajNuclc), c(1:3,4:7,10:11,8:9)]
dat[(dat$MajNuc!=dat$MinNuclc & dat$MinNuc==dat$MajNuclc),] <- dat[(dat$MajNuc!=dat$MinNuclc & dat$MinNuc==dat$MajNuclc), c(1:3,6:7,4:5,8:11) ]
dat[(dat$MajNuc!=dat$MajNuclc & dat$MinNuc==dat$MinNuclc),] <- dat[(dat$MajNuc!=dat$MajNuclc & dat$MinNuc==dat$MinNuclc), c(1:3,6:7,4:5,10:11,8:9)]
dat$GPOS <- paste0(dat$CHROM, "_", dat$POS)

# confirm same base pairs
(sb2 <- sum(dat$MajNuc == dat$MajNuclc & dat$MinNuc == dat$MinNuclc))
(sb1 <- sum(dat$MajNuc == dat$MajNuclc & dat$MinNuc != dat$MinNuclc))
(sb0 <- sum(dat$MajNuc != dat$MajNuclc & dat$MinNuc != dat$MinNuclc))
nrow(dat); (sb0 + sb1 + sb2); nrow(dat) - (sb0 + sb1 + sb2) == 0

# Write a function to rotate the .frq files from earlier.
# For each position, calculate the average major (specified) allele frequency
# and the number of alleles used in that calculation (i.e., n).
afcalc <- function(df, allele) df[,c("CHROM", "POS", "sample", allele)] %>% 
  pivot_wider(names_from = sample, values_from = allele) %>% 
  group_by(CHROM, POS) %>% 
  summarise(AF = mean(c_across(3:(ncol(.)-2)), na.rm = T))

# Run of the hdGBS dataset focusing on the major allele.
hdb_afs <- afcalc(df = dat, allele = "MajAF")
write.csv(hdb_afs, "hdgbs_afs.csv", row.names = F)
lcb_afs <- afcalc(df = dat, allele = "MajAFlc") 
write.csv(lcb_afs, "lcb_afs.csv", row.names = F)


# Together ---------------------------------------------------------------------

lc_beagle <- read.csv("lcb_afs.csv")
hd_beagle <- read.csv("hdgbs_afs.csv")
hd_mcgeno <- as.data.frame(afsm %>% 
  mutate(CHR = str_sub(GPOS, start = 0, end = 11),
         POS = as.numeric(str_sub(GPOS, start = 13))) %>% 
  select(c("CHR", "POS", "AF"))) %>% 
  rename("CHROM" = CHR)

full <- merge(lc_beagle, hd_beagle, by = c("CHROM", "POS")) %>% 
  merge(., hd_mcgeno, by = c("CHROM", "POS")) %>% 
  `colnames<-`(.,c("CHROM", "POS", "lcbeagle", "hdbeagle", "hdmcgeno")) %>% 
  pivot_longer(cols = c("hdbeagle", "hdmcgeno"))

labs <- c(`hdbeagle` = "Beagle imputation",
          `hdmcgeno` = "Most common genotype")

ggplot(data = full, 
       aes(x = lcbeagle,
           y = value)) +
  geom_smooth(method = "lm", 
              alpha = 1/6,
              colour = "black",
              linetype = 2) +
  geom_point(shape  = 21, 
             colour = "black", 
             fill   = "grey50", 
             alpha  = 1/4) + 
  facet_wrap(~name, labeller = as_labeller(labs)) + theme_bw() +
  stat_poly_eq(use_label(c("R2", "p")),
               label.x = "left",
               label.y = "top",
               small.p = TRUE) +
  labs(x = "hdGBS Allele Frequency",
       y = "Imputed lcWGS Allele Frequency")

ggsave("../../plots/hdgbs_imp_methods_lc.tiff",
       width = 12, height = 8, dpi = 300)






