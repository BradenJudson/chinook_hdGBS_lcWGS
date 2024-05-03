setwd("~/ots_landscape_genetics/lcWGS/stats")

library(tidyverse)

dat <- as.data.frame(matrix(read.table("alignment_flagstat.txt",
                     sep = "\t", col.names = "all")[,1],
                     ncol = 17, byrow = T)) %>% 
  `colnames<-`(., c("file", "reads", "primary","secondary", "supplementary",
                    "duplicates", "p.duplicates", "mapped", "p.mapped",
                    "paired", "read1", 
                    "read2", "p.paired", "self.paired", "singletons",
                    "mmap.diff", "mmap.diffq5")) %>% 
  mutate(reads = as.numeric(gsub("\\D+", "", reads))/10,
         p.paired = as.numeric(str_sub(p.paired, -13, -9)),
         file = gsub("^.*\\.","", str_sub(start = 40,end = 70, 
                                  substr(file, 1, nchar(file)-13))),
         abb = str_sub(start = 0, end = 3, file)) %>% 
  select(!c(duplicates, secondary, supplementary, mapped))

hist(dat$reads, breaks = 100, main = NULL, xlab = "Reads")

write.csv(dat[,c(1:2)], "lcwgs_indv_reads.csv", row.names = F)

sites <- read.delim(file = "../../data/ch2023_sequenced.txt") %>% 
  arrange(Latitude) %>% 
  mutate(site = tools::toTitleCase(tolower(gsub("\\_.*", "", Population))),
         sitenum = as.numeric(rownames(.)),
         abb = tools::toTitleCase(str_sub(start = 0, end = 3, site)))

fd <- merge(dat, sites, by = "abb") %>% select(!c("abb"))

ggplot(data = fd, aes(x = site, y = reads/1e6)) + 
  geom_hline(yintercept = 20, linetype = 2, colour = "red") +
  geom_boxplot(alpha = 1/10) +  theme_bw() +
  geom_point(size  = 3/2, shape = 21,
             color = "black", fill = "grey90") +
  theme(axis.text.x = element_text(angle = 90,
                      vjust = 0.5, hjust = 1),
        legend.position = "none") +
  labs(x = NULL, y = "Reads (millions)")

ggsave("../../plots/lcwgs_indv_reads.tiff", dpi = 300, height = 7, width = 12)

