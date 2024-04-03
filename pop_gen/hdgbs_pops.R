setwd("~/ots_landscape_genetics/pop_gen")

library(tidyverse); library(vcfR); library(ggplot2); library(adegenet)

hdgbs <- read.vcfR("../data/hdgbs_snps_maf2_m60.vcf.gz")

hdgl <- vcfR2genlight(hdgbs)

# hdgl_sub <- gl.drop.ind(hdgl)

# Takes a few hours.
# hd_pca <- glPca(hdgl, nf = 3, parallel = T, n.cores = 12)

# write.csv(hd_pca$scores, "hdgbs_pcascores.csv", row.names = TRUE)



"%ni%" <- Negate("%in%")
# drop_indvs <- rownames(pc_scores[pc_scores$PC1 > 50, ])
# hdgl_sub <- hdgl[indNames(hdgl) %ni% drop_indvs]
# hd_sub_pca <- glPca(hdgl_sub, nf = 3, parallel = T, n.cores = 12)
# write.csv(hd_sub_pca$scores, "hdgbs_sub_pcascores.csv", row.names = TRUE)


# Eigenvalues (% variation explained for each PC axis).
(pc_var <- c(hd_sub_pca$eig/sum(hd_sub_pca$eig)*100)[1:5]); barplot(hd_pca$eig)

coords <- read.table("../map/ch2023_sequenced.txt", sep = "\t", header = T) %>% 
  mutate(short = str_to_title(substr(Population, start = 1, stop = 3))) %>% 
  filter(Population != "TAKHANNE_RIVER")

coords[coords$Population == "TAKHANNE_RIVER", "short"] <- "Takha"
coords[coords$Population == "KILDALA_RIVER", "short"] <- "Kild"


pc_scores <- as.data.frame(hd_sub_pca$scores) %>% 
  mutate(pop = gsub("-.*","", rownames(.))) %>% 
  merge(., coords, by.x = "pop", by.y = "short") %>% 
  arrange(desc(Latitude)) %>% 
  mutate(pop = factor(reorder(pop, Latitude)))


ggplot(data = pc_scores, 
       aes(x = PC1, y = PC2,  fill = factor(Latitude), group = pop)) +
  geom_point(shape = 21, size = 3/2) +
  labs(x = paste0("PC1 (", round(pc_var[1], 1), "%)"),
       y = paste0("PC2 (", round(pc_var[2], 1), "%)")) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "right") +
  scale_fill_manual(values = alpha(c(viridis_pal(option = "D")(length(unique(pc_scores$pop)))), 3/4),
                    labels = levels(pc_scores$pop)) +
  guides(fill=guide_legend(ncol = 3, reverse = TRUE))
  

ggsave("plots/pca_fillscale.tiff", dpi = 300, width = 8, height = 5)
