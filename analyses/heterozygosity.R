setwd("~/ots_landscape_genetics/analyses/")

library(tidyverse)

pop_map <- read.csv("../data/shared_samples_n361.csv")

# Function for converting VCF to bed (usable by pcadapt).
# Assumes VCF is gzipped.
vcf2het <- \(vcf_file) system(paste("plink.exe --vcf", vcf_file, 
                                   " --het --hardy --aec --double-id --out", 
                                   gsub("\\.vcf.gz", "", vcf_file)))

# Make bed files necessary for pcadapt.
vcf2het("../data/vcfs_n361/hdgbs_maf5_m70_pruned.vcf.gz")
vcf2het("../data/vcfs_n361/lcwgs_ldpruned_maf005_n361.vcf.gz")




## NEED TO DOUBLE CHECK - WE WANT LD-PRUNED HERE? WHAT ABOUT HWE?




het <- \(het_file, id_col) {
  
  het <- read.delim(het_file, sep = "") %>% 
    mutate(nhet = N.NM.-O.HOM.,
           phet = nhet/N.NM.) %>% 
    merge(., pop_map,
          by.x = "IID",
          by.y = {{id_col}}) %>%
    group_by(site_full) %>%
    summarise(avg_ohet = mean(phet),
              sd_ohet  = sd(phet)) %>% 
    arrange(avg_ohet) %>% 
    mutate(site_full = fct_reorder(site_full, avg_ohet))
}


hdgbs_het <- het("../data/vcfs_n361/hdgbs_maf5_m70.het", id_col = "fish_ID") 
lcwgs_het <- het("../data/vcfs_n361/lcwgs_ldpruned_maf005_n361.het", id_col = "bam")

ggplot(data = hdgbs_het, aes(x = site_full, y = avg_ohet)) + 
  geom_errorbar(aes(ymin=avg_ohet-sd_ohet, ymax = avg_ohet+sd_ohet)) +
  geom_point() +
  labs(x = NULL)

hets <- merge(
  hdgbs_het %>% dplyr::rename(hd_het = avg_ohet, hd_sd = sd_ohet),
  lcwgs_het %>% dplyr::rename(lc_het = avg_ohet, lc_sd = sd_ohet),
  by = "site_full"
)

ggplot(data = hets, 
       aes(x = hd_het, y = lc_het)) +
  geom_errorbar(aes(xmin = hd_het-hd_sd,
                    xmax = hd_het+hd_sd), 
                alpha = 1/5) +
  geom_errorbar(aes(ymin = lc_het-lc_sd, 
                    ymax = lc_het+lc_sd), 
                alpha = 1/5) +
  stat_smooth(method = "lm",
              colour = "black",fullrange = T) +
  geom_point(shape = 21, size = 2, fill = "gray") +
  stat_poly_eq(use_label(c("R2")), label.x = "right",
               label.y = "bottom") +
  labs(x = "hdGBS heterozygosity estimates",
       y = "lcWGS heterozygosity estimates") +
  theme_bw() +
  theme(panel.grid = element_blank())


