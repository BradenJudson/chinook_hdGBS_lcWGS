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
# coords[coords$Population == "TAKHINI_RIVER", "short"] <- "Takhi"
coords[coords$Population == "TAKHANNE_RIVER", "short"] <- "Takha"

pc_scores <- as.data.frame(hd_sub_pca$scores) %>% 
  mutate(pop = gsub("-.*","", rownames(.))) %>% 
  merge(., coords, by.x = "pop", by.y = "short") %>% 
  arrange(desc(Latitude)) %>% 
  mutate(pop = factor(pop))

pc_scores$pop2 <- factor(pc_scores$pop, levels = sample(levels(pc_scores$Latitude)))

pc_scores$pop <- factor(pc_scores$pop, levels = pc_scores$pop[order(pc_scores$Latitude)],
                        ordered = T)

# ggplot(data = pc_scores, aes(x = PC1, y = PC2)) +
#   geom_point(aes(color = pop)) +
#   theme_bw() +
#   theme(legend.position = "none")

ggplot(data = pc_scores, 
       aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = reorder(pop, Latitude)), shape = 21, size = 2.5) +
  labs(x = paste0("PC1 (", round(pc_var[1], 1), "%)"),
       y = paste0("PC2 (", round(pc_var[2], 1), "%)")) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "top") +
  guides(fill = guide_legend(nrow = 5)) +
  scale_fill_manual(values = c("#f0f921", "#f2f227", "#f5eb27", "#f7e225",
                               "#f9dc24", "#fbd324", "#fccd25", "#fdc627",
                               "#febe2a", "#feb82c", "#fdb130", "#fdab33",
                               "#fca537", "#fa9e3b", "#f9983e", "#f79143",
                               "#f58c46", "#f3854b", "#f0804e", "#ee7b51",
                               "#eb7556", "#e87059", "#e56a5d", "#e26561",
                               "#de6164", "#da5b69", "#d7566c", "#d35171",
                               "#cf4c74", "#cc4778", "#c7427c", "#c33d80",
                               "#be3885", "#ba3388", "#b42e8d", "#b02991",
                               "#ab2494", "#a51f99", "#a01a9c", "#99159f",
                               "#9410a2", "#8e0ca4", "#8707a6", "#8104a7",
                               "#7a02a8", "#7401a8", "#6c00a8", "#6600a7",
                               "#6001a6", "#5801a4", "#5102a3", "#4903a0",
                               "#43039e", "#3c049b", "#330597", "#2c0594",
                               "#220690", "#19068c", "#0d0887"))

  
  
  
ggsave("plots/pca_fillscale.tiff", dpi = 300, width = 10, height = 7)