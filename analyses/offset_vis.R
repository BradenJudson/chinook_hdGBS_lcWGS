setwd("~/ots_landscape_genetics/analyses")

library(tidyverse); library(sf); library(vcfR); library(ggpmisc)
library(cowplot); library(viridis); library(ggpubr); library(ConGenFunctions)

# Data for each location and coordinates.
dat <- read.csv("../data/chinook_full_sample_info.csv")
loc <- read.csv("../data/ch_site_info.csv")[,c("Site", "Latitude", "Longitude")]

sites <- merge(dat, loc, by.x = "site_full", by.y = "Site")

# Heterozygosity and offset estimates for each population.
het <- read_delim("../data/offset/het_estimates.txt")
off85 <- read.delim("../data/offset/genomic_offsets.txt")[,c(1,3)]
off26 <- read.delim("../data/offset/genomic_offsets.txt")[,c(1:2)]

# Join site info with offset estimates.
data <- merge(sites, off85, by.x = "site_full", by.y = "population") %>% 
  merge(., off26, by.x = "site_full", by.y = "population")

# Samples per population.
sample_counts <- sites %>% 
  group_by(site_full) %>% 
  tally()

summary(sample_counts$n)

# Combine site info and heterozygosity.
# Keep one value per population.
hetoff <- cbind(het, off85)[,2:7] %>% 
  merge(., sites[,c(1,4:5)], 
        by.x = "population", 
        by.y = "site_full") %>% 
  group_by(population) %>% 
  dplyr::sample_n(1)

(hetlm <- ggplot(data = hetoff,
       aes(x = indv_het_mean,
           y = offset85)) +
  stat_smooth(method = "lm",
              alpha  = 1/6,
              colour = "black",
              linewidth = 1,
              lineend = "round") +
  stat_poly_eq(use_label(c("R2", "p")),
               label.x = "left",
               label.y = "top",
               small.p = TRUE) + theme_bw() + 
  geom_point(size = 2, shape = 21, aes(fill = as.factor(Latitude))) +
  scale_fill_manual(values = c(viridis_pal(option = "D")(length(unique(hetoff$population)))),
                    labels = hetoff$population)  +
    labs(y = "Genomic offset",
         x = "Mean individual heterozygosity") +
  guides(fill = guide_legend(ncol = 2, byrow = TRUE)) +
  theme(legend.title = element_blank()))

ggsave("../plots/offset_het.tiff", dpi = 300,
       width = 11, height = 8, bg = 'white')
  

# Visualize w/ map -------------------------------------------------------------


USA <- sf::st_as_sf(geodata::gadm(country = "USA", level = 0, path = "../map/"))
CAN <- sf::st_as_sf(geodata::gadm(country = "CAN", level = 1, path = "../map/"))

hets <- cbind(loc[!loc$Site %in% c("Raft", "Harrison", "Nahatlatch"),], het)

offset_plot <- \(df, variable, midpoint, cgr_rev, plot_title) {
  
  p1 <- ggplot(data = df) + theme_bw() + 
    geom_sf(data = USA, fill = "gray90", linewidth = 1/10) +
    geom_sf(data = CAN, fill = "gray90", linewidth = 1/10) +
    coord_sf(xlim = c(-115, -165), ylim = c(41, 66)) +
    labs(x = NULL, y = NULL) +
    geom_point(aes(x = Longitude, 
                   y = Latitude, 
                   fill = {{variable}}),
               shape = 21, stroke = 1/10) +
    theme(legend.background = element_rect(colour = 'black', fill = 'white'),
          legend.position = "right", legend.title = element_blank()) +
    scale_fill_gradient2(low  = "skyblue", 
                         high = "red2", 
                         midpoint = midpoint,
                         guide = guide_colorbar(draw.ulim = TRUE, 
                                                draw.llim = TRUE, 
                                                reverse = cgr_rev,
                                                ticks   = TRUE,
                                                frame.colour = 'black',
                                                ticks.color  = 'black')) +
    ggtitle(plot_title)
  
  p1_leg <- as_ggplot(cowplot::get_legend(p1))
  
  p2 <- p1 + theme(legend.position = "none") +
      geom_point(aes(x = Longitude, 
                     y = Latitude, 
                     fill = {{variable}},
                     size = {{variable}}),
                 shape = 21)


  ConGenFunctions::insettr(p2, p1_leg,
                           location = "tr",
                           height = 0.3,
                           width = 0.1)
  
  
  
}

(off85 <- offset_plot(data, variable = offset85, 
                      midpoint = mean(data$offset85), 
                      cgr_rev = FALSE,
                      plot_title = "Genomic offset (ssp85)") )
ggsave("../plots/genomic_offset85.tiff", bg = 'white',
       dpi = 300, width = 8, height = 6)

(off26 <- offset_plot(data, variable = offset26, 
                      midpoint = 0.059, 
                      cgr_rev = FALSE,
                      plot_title = "Genomic offset (ssp26)"))
ggsave("../plots/genomic_offset26.tiff", bg = 'white',
       dpi = 300, width = 8, height = 8)

# https://ggplot2.tidyverse.org/reference/guide_colourbar.html
(hetpl <- offset_plot(hets, 
                      variable = indv_het_mean, 
                      midpoint = 0.24, 
                      cgr_rev = TRUE,
                      plot_title = "Average individual heterozygosity"))
ggsave("../plots/indv_heterozygosity.tiff", bg = 'white',
       dpi = 300, width = 8, height = 8)

cowplot::plot_grid(plotlist = list(off85, hetpl), ncol = 2)
ggsave("../plots/offset_het.tiff", width = 14, height = 8, dpi = 300)


