setwd("~/ots_landscape_genetics/comp")

library(tidyverse); library(ggpmisc)

hdgbs_reads <- read.table("../hdGBS/stats/unfiltered_stats/indv_reads.txt",
                          sep = "") %>% 
  `colnames<-`(., c("Sample", "hdgbs_reads")) %>% 
  mutate(Sample = gsub("\\.", "-",gsub("\\.1.fq.gz", "", Sample)))
  

lcwgs_reads <- read.csv("../lcWGS/stats/lcwgs_indv_reads.csv") %>% 
  `colnames<-`(., c("Sample", "lcwgs_reads")) 
colnames(hdgbs_reads); colnames(lcwgs_reads)

com <- merge(hdgbs_reads, lcwgs_reads, by = "Sample")

lcwgs_reads[!lcwgs_reads$Sample %in% c(hdgbs_reads$Sample),]
hdgbs_reads[!hdgbs_reads$Sample %in% c(lcwgs_reads$Sample),]

ggplot(data = com, 
       aes(x = hdgbs_reads/1e6, 
           y = lcwgs_reads/1e6)) + 
  geom_smooth(method = "lm",
              colour = "red",
              alpha = 1/6,
              linetype= "dashed") +
  geom_point(alpha = 4/5, shape = 21,
             fill = "gray80", size = 2,
             colour = "black") +
  theme_bw() +
  labs(x = "hdGBS Reads (millions)",
       y = "lcWGS Reads (millions)") +
  stat_poly_eq(use_label(c("R2", "p")),
               label.x = "right",
               label.y = "top",
               small.p = T)



# Ridgeline plots --------------------------------------------------------------

# sample_info <- read.csv("../hdGBS/info_files/hdGBS_sampleinfo.csv") 

lf <- com %>% pivot_longer(cols = c("lcwgs_reads", "hdgbs_reads"),
                           names_to = "seq", values_to = "reads") %>% 
  mutate(pop = factor(str_sub(start = 0, end = 3, Sample)),
         seq = as.factor(seq),
         lab = case_when(seq == "hdgbs_reads" ~ "hdGBS",
                         seq == "lcwgs_reads" ~ "lcWGS"))

(gr <- ggplot(data = lf, 
       aes(y = pop, 
           x = reads/1e6,
           fill = seq,
           colour = seq)) +
  geom_density_ridges(alpha = 3/5,
                      scale = 9/10,
                      jittered_points = T, 
                      point_shape = "|", point_size = 3/2,
                      position = position_points_jitter(height = 0),
                      rel_min_height = 0.01) +
  scale_fill_manual(values = c("#D55E0050", "#0072B250"), 
                    labels = c("hdGBS", "lcWGS")) +
  scale_color_manual(values = c("#D55E00", "#0072B2"), guide = "none") +
  theme_bw() +
  labs(x = "Reads (millions)",
       y = NULL) +
  scale_x_continuous(breaks = seq(0, 150, 25)) +
  theme(legend.position = "top",
        legend.title = element_blank()))

ggsave("../plots/read_dists.tiff", dpi = 300,
       height = 16, width = 9)


# Boxplots ---------------------------------------------------------------------


(bp <- ggplot(data = lf, 
       aes(x = pop, 
           y = reads/1e6)) +
  geom_boxplot(alpha = 1/10) +
  facet_wrap(~lab, ncol = 1) +
  geom_point(size = 3/2,
             shape = 21,
             color = "black",
             fill = "gray80") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90,
                                   vjust = 0.5, 
                                   hjust = 1)) +
  labs(x = NULL, y = "Reads (millions)"))

ggsave("../plots/read_boxplots.tiff", dpi = 300,
       height = 10, width = 14)

(hi <- ggplot(data = lf, 
              aes(x = reads/1e6)) +
    geom_histogram(colour = "black",
                   fill = "gray80",
                   bins = 100) +
    theme_bw() +
    facet_wrap(~lab, ncol = 1) +
    scale_y_continuous(position = "right") +
    labs(x = "Reads (millions)", y = "Count") +
    scale_x_continuous(breaks = seq(0, 120, 20))) 


cowplot::plot_grid(bp, hi, align = "VH", ncol = 2, rel_widths = c(1, 0.8))

ggsave("../plots/read_plots.tiff", dpi = 300,
       height = 12, width = 16)




