setwd("~/ots_landscape_genetics/lcWGS/bam_pop_sep")

library(tidyverse)

bams <- read.csv("../../data/shared_samples_n361.csv") %>% 
  group_by(site_full) %>% 
  split(., f = .$site_full) %>% 
  lapply(., \(x) x[,1])

mapply(write.table, bams, file = paste0("n361_bams/", tolower(names(bams)), "_bams.txt"),
       MoreArgs = list(quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\n"))

