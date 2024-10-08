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

# Shapefiles for land.
USA <- sf::st_as_sf(geodata::gadm(country = "USA", level = 0, path = "../map/"))
CAN <- sf::st_as_sf(geodata::gadm(country = "CAN", level = 1, path = "../map/"))

hets <- cbind(loc[!loc$Site %in% c("Raft", "Harrison", "Nahatlatch"),], het)

# Function for plotting offset/diversity metrics.
offset_plot <- \(df, variable, midpoint, cgr_rev, plot_title) {
  
  # For consistent label bins.
  qs <- round(as.numeric(quantile(df[,deparse(substitute(variable))])), 2)
  
  # Plot land and populations on it.
  p1 <- ggplot(data = df) + theme_bw() +                                    
    geom_sf(data = USA, fill = "gray90", linewidth = 1/10) +
    geom_sf(data = CAN, fill = "gray90", linewidth = 1/10) +
    coord_sf(xlim = c(-115, -165), ylim = c(41, 66)) +
    labs(x = NULL, y = NULL) +
    geom_point(aes(x = Longitude, 
                   y = Latitude, 
                   fill = {{variable}},   # Colour is metric of interest.
                   size = {{variable}}),  # Point size is also proportional to value.
               shape = 21, stroke = 1/10) + ggtitle(plot_title) + 
    theme(legend.background = element_rect(colour = 'black', fill = 'white'),
          legend.position = "right", legend.title = element_blank(),
          legend.justification = "top") +
    scale_fill_gradient2(low  = if(cgr_rev == FALSE) "skyblue" else "red2", 
                         high = if(cgr_rev == FALSE) "red2" else "skyblue", 
                         midpoint = qs[[3]], breaks = qs, 
                         limits = c(min(qs), max(qs))) +
    scale_size_continuous(range = c(1,5), breaks = c(qs),
                          limits = c(min(qs), max(qs))) +
    guides(fill = guide_legend(), size = guide_legend()) 
    # Above ensures size and fill legends are combined.

   # Extract legend from plot and remove margins.  
   p1_leg <- as_ggplot(cowplot::get_legend(p1)) +
    theme(plot.margin = unit(c(0,0,0,0), "points"))
  
  # Remove legend from main plot.
  p2 <- p1 + theme(legend.position = "none")
  
  # Use custom function to inset legend inside figure.
  ConGenFunctions::insettr(p2, p1_leg,
                           location = "tr",
                           height = 0.3,
                           width = 0.1)
   
}

(off85 <- offset_plot(data, variable = offset85, 
                      midpoint = mean(data$offset85), 
                      cgr_rev = FALSE,
                      plot_title = "Genomic offset (ssp85)"))
ggsave("../plots/genomic_offset85.tiff", bg = 'white',
       dpi = 300, width = 8, height = 6)

(off26 <- offset_plot(data, variable = offset26, 
                      midpoint = 0.059, 
                      cgr_rev = FALSE,
                      plot_title = "Genomic offset (ssp26)"))
ggsave("../plots/genomic_offset26.tiff", bg = 'white',
       dpi = 300, width = 8, height = 6)

(hetpl <- offset_plot(hets, 
                      variable = indv_het_mean, 
                      midpoint = 0.24, 
                      cgr_rev = TRUE,
                      plot_title = "Average individual heterozygosity"))
ggsave("../plots/indv_heterozygosity.tiff", bg = 'white',
       dpi = 300, width = 8, height = 6)

cowplot::plot_grid(plotlist = list(off85, hetpl), ncol = 2)
ggsave("../plots/offset_het.tiff", width = 14, height = 8, dpi = 300)


