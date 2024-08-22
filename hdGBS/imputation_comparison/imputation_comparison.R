setwd("~/ots_landscape_genetics/hdGBS/imputation_comparison")

library(vcfR); library(tidyverse)

# Part I: Replace missing with most common genotype ----------------------------

# Read in original VCF with missing data.
original <- read.vcfR("vcfs/hdgbs_subset_134kSNPs_n362_original.vcf")
o.mat <- extract.gt(original) # Extract genotypes.

sum(is.na(o.mat)) # Number of missing sites.
round(sum(is.na(o.mat))/(nrow(o.mat)*ncol(o.mat))*100, 2) # Proportion of missing sites.

missing <- which(is.na(o.mat), T) # Find positions of all missing SNPs.

# Replace missing genotypes with most common genotype at that site.
mcg <- apply(o.mat, 2, function(x) replace(x, is.na(x), names(which.max(table(x))))) 

# Isolate imputed genotypes using this approach.
mcsub <- mcg[missing]

# Part II: Use BEAGLE to imputed missing genotypes -----------------------------

# Read in VCF where missing data were imputed using BEAGLE v5.
imputed <- read.vcfR("vcfs/hdgbs_subset_134kSNPs_n362_imputed.vcf")
i.mat <- extract.gt(imputed) # Extract genotypes.

# Isolate imputed genotypes only.
isub <- i.mat[missing]; length(isub)

# Confirm ordering and dimensions are consistent.
all.equal(dim(mcg), dim(o.mat), dim(i.mat))
all.equal(length(missing)/2, length(mcsub), length(isub))

# Combine everything into a single dataframe.
# Requires some reformatting.
full <- as.data.frame(missing) %>% 
  mutate(beagle = sub("\\|", "/", isub),
         mostco = mcsub)

# Variation in imputation methods.
(match  <- sum(full$beagle == full$mostco))  # 8,222,927
(nmatch <- sum(full$beagle != full$mostco))  # 1,004,461
(total  <- nrow(full))                       # 9,227,388
(prop   <- round(match/total*100, 2))        # 89.11%
