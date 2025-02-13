setwd("~/ots_landscape_genetics/analyses")

library(tidyverse); library(cowplot); library(ggpmisc); library(ggExtra)


# lcWGS data -------------------------------------------------------------------

lcwgs_reads <- read_delim("../lcWGS/stats/reads.txt", col_names = c("sample", "reads"))
lcwgs_depth <- read_delim("../lcWGS/stats/bam_coverage.txt", col_names = c("sample", "depth"))

lcwgs <- merge(lcwgs_depth, lcwgs_reads, by = "sample")

n361 <- read.delim("../data/pop_map_n361.txt", col.names = c("sample", "pop"))


# hdGBS data -------------------------------------------------------------------

hdgbs_depth <- read_delim("../hdGBS/stats/unfiltered_stats_complete/out.idepth") %>% 
  rename("Sample" = INDV)
hdgbs_reads <- read.csv("../hdGBS/stats/sample_reads.csv", row.names = 1)[,c(1:2)]

hdgbs <- merge(hdgbs_depth, hdgbs_reads) %>% filter(!Sample %in% n361$sample)

# ------------------------------------------------------------------------------

vis <- \(df, x, y) {
  p1 <- ggplot(data = df,
         aes(x = {{x}}/1e6,
             y = {{y}})) +
    geom_segment(aes(y = mean({{y}}), 
                     yend = mean({{y}}), 
                     x = -Inf, xend = Inf),
                 linewidth = 1/2, linetype = 2) +
    geom_segment(aes(x = mean({{x}})/1e6,
                     xend = mean({{x}})/1e6,
                     y = -Inf, yend = Inf),
                 linewidth = 1/2, linetype = 2) +
    geom_smooth(method = "lm",
                colour = "red",
                alpha  = 1/6,
                linetype = 1) +
    geom_point(colour = "black",
               shape  = 21, 
               fill   = "gray80",
               size   = 2,
               alpha  = 4/5) +
    theme_bw() +
    theme(panel.grid = element_line(color = "gray95"),
          panel.grid.minor.y = element_blank()) +
    labs(x = "Reads (millions)", 
         y = "Average Individual Coverage") + 
    stat_poly_eq(use_label(c("R2", "p")),
                 label.x = "left",
                 label.y = "top",
                 small.p = TRUE)
  
  (p2 <- ggMarginal(p1, type = "density",
                    colour = "black",
                    fill = "grey90",
                    alpha = 7/10,
                    linewidth = 1))
    
}

(hd <- vis(hdgbs, x = Reads, y = MEAN_DEPTH))
(lc <- vis(df = lcwgs, x = reads, y = depth))

cowplot::plot_grid(plotlist = list(hd, lc),
                   ncol = 1, labels = c("a)", "b)"))

ggsave("../plots/reads_depth.tiff", dpi = 300, 
       bg = 'white', width = 12, height = 12)

