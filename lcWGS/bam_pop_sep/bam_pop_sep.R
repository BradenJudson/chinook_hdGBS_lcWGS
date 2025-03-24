setwd("~/ots_landscape_genetics/lcWGS/bam_pop_sep")

library(tidyverse)

bams <- read.csv("../../data/n385_shared_indvs.csv", 
                 col.names = c("id", "bam", "pop"),
                 header = FALSE) %>% 
  group_by(pop) %>% 
  split(., f = .$pop) %>% 
  lapply(., \(x) x[,"bam"])

mapply(write.table, bams, file = paste0("n385_bams/", tolower(names(bams)), "_bams.txt"),
       MoreArgs = list(quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\n"))

