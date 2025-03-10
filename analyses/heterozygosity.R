setwd("~/ots_landscape_genetics/analyses/")

library(tidyverse); library(ggpmisc)

pop_map <- read.csv("../data/shared_samples_n361.csv")

# Function for converting VCF to bed (usable by pcadapt).
# Assumes VCF is gzipped.
vcf2het <- \(vcf_file) system(paste("plink.exe --vcf", vcf_file, 
                                   " --het --aec --double-id --out", 
                                   gsub("\\.vcf.gz", "", vcf_file)))

# Make bed files necessary for pcadapt.
vcf2het("../data/vcfs_n361/hdgbs_maf5_m70.vcf.gz")
vcf2het("../data/vcfs_n361/lcwgs_ldpruned_maf005_n361.vcf.gz")
vcf2het("../data/vcfs_n361/chinook_filtered_maf5_imputed_commonsubset.vcf.gz")
vcf2het("../data/vcfs_n361/hdgbs_n361_maf5_m70_commonsubset.vcf.gz")

# https://www.biostars.org/p/266502/

het <- \(het_file, id_col) {
  
  # observed / non-missing * 100 (TH)
  
  het <- read.delim(het_file, sep = "") %>% 
    mutate(nhet = N.NM.-O.HOM.,
           phet = nhet/N.NM.) %>% 
    merge(., pop_map,
          by.x = "IID",
          by.y = {{id_col}}) %>%
    group_by(IID) %>%
    summarise(avg_ohet = mean(phet),
              sd_ohet  = sd(phet)) %>% 
    arrange(desc(avg_ohet))
    # mutate(site_full = fct_reorder(site_full, avg_ohet),
    #        rank = as.numeric(rownames(.)))
}


hdgbs_het <- het("../data/vcfs_n361/hdgbs_maf5_m70.het", id_col = "fish_ID") 
lcwgs_het <- het("../data/vcfs_n361/lcwgs_ldpruned_maf005_n361.het", id_col = "bam")

lcwgs_sub <- het("../data/vcfs_n361/chinook_filtered_maf5_imputed_commonsubset.het", id_col = "bam") %>% 
  merge(., pop_map, by.x = "IID", by.y = "bam")

hdgbs_sub <- het("../data/vcfs_n361/hdgbs_n361_maf5_m70_commonsubset.het", id_col = "fish_ID")

ggplot(data = hdgbs_het, aes(x = site_full, y = avg_ohet)) + 
  geom_errorbar(aes(ymin=avg_ohet-sd_ohet, ymax = avg_ohet+sd_ohet)) +
  geom_point() +
  labs(x = NULL)

hets <- merge(
  hdgbs_sub %>% dplyr::rename(hd_het = avg_ohet, hd_sd = sd_ohet),
  lcwgs_sub %>% dplyr::rename(lc_het = avg_ohet, lc_sd = sd_ohet),
  by.x = "IID", by.y = "fish_ID"
)

(reg <- ggplot(data = hets, 
       aes(x = hd_het, y = lc_het)) +
  # geom_errorbar(aes(xmin = hd_het-hd_sd,
  #                   xmax = hd_het+hd_sd), 
  #               alpha = 1/5) +
  # geom_errorbar(aes(ymin = lc_het-lc_sd, 
  #                   ymax = lc_het+lc_sd), 
  #               alpha = 1/5) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
  stat_smooth(method = "lm",
              colour = "black",fullrange = T) +
  geom_point(shape = 21, size = 2, fill = "gray") +
  stat_poly_eq(use_label(c("R2")), label.x = "right",
               label.y = "bottom") +
  labs(x = "hdGBS ranked heterozygosity estimates",
       y = "Imputed ranked lcWGS heterozygosity estimates") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_y_continuous(transform = "reverse") +
  scale_x_continuous(transform = "reverse") )

ggplotly(reg)


# test --------------------------------------------------------------------

library(vcfR)

lc <- read.vcfR("../data/vcfs_n361/hdgbs_n361_maf5_m70_commonsubset.vcf.gz")

