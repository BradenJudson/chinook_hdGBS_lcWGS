setwd("~/ots_landscape_genetics/pop_gen")

library(tidyverse); library(ggplot2)

covmat <- read.table("../data/pca_m15_maf005.cov")
mmepca <- eigen(covmat)
eigvec <- mmepca$vectors %>% as.data.frame()

files <- read.table("../data/lcwgs_bam_list_n453.txt", col.names = c("file")) %>% 
  mutate(fish_ID = gsub(".dedup.clip.bam", "", 
                        gsub("^.*IDT_i\\d{1,3}_\\d{1,3}\\.", "", file))) 

sites <- read.delim("../data/ch2023_sequenced.txt")[,c(4,5,8:9L)] %>% 
  rename("site_full" = "Site")
sites[sites$site_full == "Pitt", "site_full"] <- "Upper Pitt"

samples <- read.csv("../data/landgen_chinook_indvs.csv")[,c(1,3)] %>% 
  group_by(fish_ID) %>% sample_n(1) %>% 
  mutate(fish_ID = gsub("\\.", "-", fish_ID),
         site_full = gsub("[[:space:]]Brood", "", site_full))

pca.eigenval.sum <- sum(mmepca$values)
varPC1 <-  (mmepca$values[1]/pca.eigenval.sum)*100
varPC2 <-  (mmepca$values[2]/pca.eigenval.sum)*100

pcavec <- as_tibble(cbind(files, eigvec)) %>% 
  group_by(fish_ID) %>% 
  merge(., samples, by = "fish_ID") %>% 
  select(c(1:5, ncol(.)))

# Get population-level average PC values for visualization.
pc_popav <- pcavec %>% 
  group_by(site_full) %>% 
  summarise(V1A = mean(V1),
            V2A = mean(V2))

# Add pop-averages to main dataframe.
pcavectors2 <- merge(pcavec, pc_popav, by = "site_full", all.x = T) %>% 
  merge(., sites, by = "site_full")

listlats <- list(North = max, South = min)
listvals <- lapply(listlats, function(x) x(pcavectors2$Latitude))

(gpca <- ggplot(data = pcavectors2, aes(fill = Latitude, colour = Latitude)) +
  geom_segment(aes(x = V1A, y = V2A, xend = V1, yend = V2), linewidth = 3/4) +
  geom_point(aes(x = V1, y = V2),color = "black", shape = 21, size = 2) + theme_bw() +
  scale_color_gradient(guide = 'none') +
  scale_fill_gradient(breaks = unlist(listvals)) +
  ggrepel::geom_label_repel(data = pc_popav, max.overlaps = Inf, 
                            min.segment.length = 0, box.padding = 1/2,
                            aes(x = V1A, y = V2A, label = site_full), 
                            inherit.aes = FALSE) +
  labs(x = paste0("PC1 (", round(varPC1, 1), "%)"),
       y = paste0("PC2 (", round(varPC2, 1), "%)")) + 
  theme(legend.position = "top", legend.title = element_blank(),
        legend.text = element_text(size = 11)) +
  guides(fill = guide_colorbar(barwidth = 25, label = T, reverse = T, 
                               ticks = F, frame.colour = "black", label.position = "top")))


outliers <- pcavectors2[pcavectors2$V1 > 0.05, c(1,2,4)]

ggsave("../plots/lcwgs_pca.tiff", dpi = 300, width = 12, height = 9)



# Exploration ------------------------------------------------------------------

pc_var <- pcavec %>% 
  group_by(site_full) %>% 
  summarise(V1A = var(V1),
            V2A = var(V2)) %>% 
  mutate(lab = case_when(V2A > 6e-4 ~ site_full,
                         V1A > 1e-5 ~ site_full))

ggplot(data = pc_var,
       aes(x = V1A, y = V2A)) +
  geom_point() + theme_bw() +
  ggrepel::geom_label_repel(aes(label = lab),
                            max.overlaps = Inf, 
                            min.segment.length = 0) +
  labs(x = paste0("PC1 (", round(varPC1, 1), "%) Variance"),
       y = paste0("PC2 (", round(varPC2, 1), "%) Variance")) + 
  theme(legend.position = "top", legend.title = element_blank(),
        legend.text = element_text(size = 11)) 

ggsave("../plots/lcwgs_pca_variance.tiff", dpi = 300, width = 12, height = 9)



