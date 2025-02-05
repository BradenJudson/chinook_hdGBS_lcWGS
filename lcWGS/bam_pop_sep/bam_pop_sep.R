setwd("~/ots_landscape_genetics/lcWGS/bam_pop_sep")

library(tidyverse)

bams <- read.csv("../../data/shared_samples_n361.csv") %>% 
  group_by(site_full) %>% 
  split(., f = .$site_full) %>% 
  lapply(., \(x) x[,1])

mapply(write.table, bams, file = paste0("n361_bams/", tolower(names(bams)), "_bams.txt"),
       MoreArgs = list(quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\n"))


# 
# # Writes a text file containing bam files per population per line.
# # First chunk extracts sample identity.
# files <- read.table("../../data/bam_list_n453.txt", col.names = c("file")) %>% 
#   mutate(fish_ID = gsub(".dedup.clip.bam", "", 
#                         gsub("^.*IDT_i\\d{1,3}_\\d{1,3}\\.", "", file))) 
# 
# samples <- read.csv("../../data/landgen_chinook_indvs.csv") %>% 
#   mutate(fish_ID = gsub("\\.", "-", fish_ID)) %>% 
#   filter(fish_ID %in% files$fish_ID)
# 
# dat <- merge(files, samples, by = "fish_ID") %>% 
#   mutate(siteName  = tolower(gsub(" ", "_", site_full)),
#          direcfile = file) %>% 
#   group_by(site_full) %>% 
#   split(., f = .$siteName) %>% 
#   lapply(., function(x) x[,"direcfile"])
# 
# (names(dat) <- gsub("_brood", "", names(dat)))
# 
# mapply(write.table, dat, file=paste0(tolower(names(dat)), "_bams.txt"),
#        MoreArgs = list(quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\n"))

print(c(names(dat)), quote = F) # Check names.
