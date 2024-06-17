setwd("~/ots_landscape_genetics/comp")

library(tidyverse)

files <- read.table("../data/bam_list_n453.txt", col.names = c("file")) %>%
  mutate(fish_ID = gsub(".dedup.clip.bam", "",
                        gsub("^.*IDT_i\\d{1,3}_\\d{1,3}\\.", "", file))) 

hdgl <- read.table("./data/original_gls.txt",
                   col.names = c("chrom", "pos", c(files[,2]))) %>% 
  pivot_longer(cols = c(3:ncol(.))) %>% 
  separate(value, c("AA", "AB", "BB"), sep = ",") %>% 
  mutate_at(c((ncol(.)-2):ncol(.)), as.numeric)





