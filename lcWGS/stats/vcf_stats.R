setwd("~/ots_landscape_genetics/lcWGS/stats")

library(tidyverse); library(ggpmisc)


# SNP Depth --------------------------------------------------------------------


snp_depth <- read.table("../../data/lcwgs_maf5_snpdepth.mean", header = TRUE)
snp_miss <- read.table("../../data/vcf_maf005.lmiss", header = T)

snp_data <- merge(snp_depth[,c(1:3)],
                  snp_miss[,c(1:2,6)],
                  by = c(1:2)) %>% 
  pivot_longer(cols = c("MEAN_DEPTH", "F_MISS"))

rm(snp_depth); rm(snp_miss)

grid_labs <- c(`F_MISS` = "SNP Missingness (%)",
               `MEAN_DEPTH` = "SNP Depth")

ggplot(data = snp_data, aes(x = value)) +
  geom_histogram(colour = "black",
                 fill = "gray80",
                 bins = 70) +
  facet_wrap(~ name, ncol = 2, 
             scales = "free_x",
             labeller = as_labeller(grid_labs)) +
  labs(x = NULL, y = "Frequency") +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = NA,
                                    linetype = 0),
    strip.text = element_text(size = 10)
  ) 

ggsave("../../plots/lcwgs_snp_summ.png", width = 3300, height = 2000, units = "px")


(snp_summ <- snp_data %>% 
    group_by(name) %>% 
    summarise(mean = mean(value),
              med  = median(value),
              sd   = sd(value),
              min  = min(value),
              max  = max(value)) %>% 
    mutate_if(is.numeric, round, 3))


