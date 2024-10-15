setwd("~/ots_landscape_genetics/analyses")

library(tidyverse); library(sf); library(vcfR); library(ggpmisc); library(ggpmisc)
library(cowplot); library(viridis); library(ggpubr); library(ConGenFunctions)

# Data for each location and coordinates.
dat <- read.csv("../data/chinook_full_sample_info.csv")
loc <- read.csv("../data/ch_site_info.csv")[,c("Site", "Latitude", "Longitude")]

sites <- merge(dat, loc, by.x = "site_full", by.y = "Site")

# Heterozygosity and offset estimates for each population.
het <- read_delim("../data/offset/het_estimates.txt")
off85 <- read.delim("../data/offset/genomic_offsets.txt")[,c(1,3)]
off26 <- read.delim("../data/offset/genomic_offsets.txt")[,c(1:2)]
all.equal(nrow(het), nrow(off85), nrow(off26))

data <- merge(off85, off26) %>% cbind(., het[,2:ncol(het)]) %>% 
  merge(loc[loc$Site %in% off85$population,], by.x = "population", by.y = "Site")

# Combine site info and heterozygosity.
# Keep one value per population.
hetoff <- cbind(het, off85)[,2:7] %>% 
  merge(., sites[,c(1,4:5)], 
        by.x = "population", 
        by.y = "site_full") %>% 
  group_by(population) %>% 
  dplyr::sample_n(1)

(hetlm <- ggplot(data = data,
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

ggsave("../plots/offset_het_lm.tiff", dpi = 300,
       width = 11, height = 8, bg = 'white')

# Visualize w/ map -------------------------------------------------------------

# Shapefiles for land.
USA <- sf::st_as_sf(geodata::gadm(country = "USA", level = 0, path = "../map/"))
CAN <- sf::st_as_sf(geodata::gadm(country = "CAN", level = 1, path = "../map/"))

# Function for plotting offset/diversity metrics.
# Function for plotting offset/diversity metrics.
offset_plot <- \(df, variable, cgr_rev, plot_title, abs_cols, brks) {
  
  # Same as above, but for a shared/independent middle value for the legend.
  if(abs_cols == TRUE) mp <- median(c(df$offset26, df$offset85)) else
    mp <- mean(as.numeric(df[,deparse(substitute(variable))])) 
  
  # Plot land and populations on it.
  p1 <- ggplot(data = df) + theme_bw() +                                    
    geom_sf(data = USA, fill = "gray90", linewidth = 1/10) +
    geom_sf(data = CAN, fill = "gray90", linewidth = 1/10) +
    coord_sf(xlim = c(-115, -165), ylim = c(41, 66)) +
    labs(x = NULL, y = NULL) +
    geom_point(aes(x = Longitude, 
                   y = Latitude, 
                   fill  = {{variable}},   # Colour is metric of interest.
                   size  = {{variable}}),  # Point size is also proportional to value.
               shape = 21, stroke = 1/2) + ggtitle(plot_title) + 
    theme(legend.background = element_rect(colour = 'black', fill = 'white'),
          legend.position = "right", 
          legend.title = element_blank(),
          legend.justification = "top", 
          panel.grid = element_line(linetype = 3)) +
    scale_fill_gradient2(low  = if(cgr_rev == FALSE) "skyblue" else "red2", 
                         high = if(cgr_rev == FALSE) "red2" else "skyblue", 
                         midpoint = mp, 
                         breaks   = scales::breaks_pretty(n = brks),
                         limits   = if(abs_cols == TRUE) c(min(df[,"offset26"]), max(df[,"offset85"])) else
                           c(min(df[,deparse(substitute(variable))]),
                             plyr::round_any(max(df[,deparse(substitute(variable))]), 
                                             accuracy = 0.01, f = ceiling))) +
    scale_size_continuous(breaks = scales::breaks_pretty(n = brks),
                          limits = if(abs_cols == TRUE) c(min(df[,"offset26"]), max(df[,"offset85"])) else
                            c(min(df[,deparse(substitute(variable))]),
                              plyr::round_any(max(df[,deparse(substitute(variable))]), 
                                              accuracy = 0.01, f = ceiling))) +
    guides(fill = guide_legend(), 
           size = guide_legend()) 
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

# Plot worst-case offset values.
(off85 <- offset_plot(data, variable = offset85, 
                      cgr_rev = F, abs_cols = F, brks = 5,
                      plot_title = "Genomic offset (ssp85)"))
ggsave("../plots/genomic_offset85.tiff", bg = 'white',
       dpi = 300, width = 8, height = 6)

# Plot best-case offset values.
(off26 <- offset_plot(data, variable = offset26, 
                      cgr_rev = FALSE, abs_cols = F, brks = 5,
                      plot_title = "Genomic offset (ssp26)"))
ggsave("../plots/genomic_offset26.tiff", bg = 'white',
       dpi = 300, width = 8, height = 6)

# Heterozygosity plot.
(hetpl <- offset_plot(data, 
                      variable = indv_het_mean, 
                      abs_cols = F, brks = 6,
                      cgr_rev = TRUE,
                      plot_title = "Average individual heterozygosity"))
ggsave("../plots/indv_heterozygosity.tiff", bg = 'white',
       dpi = 300, width = 8, height = 6)

# Plot worst offset and heterozygosities beside each other.
cowplot::plot_grid(plotlist = list(off85, hetpl), ncol = 2)
ggsave("../plots/offset_het.tiff", width = 14, height = 8, dpi = 300)

# Plot offset values for best/worst case scenarios.
cowplot::plot_grid(plotlist = list(offset_plot(data, variable = offset26, 
                                               cgr_rev = FALSE, abs_cols = T, brks = 6,
                                               plot_title = "Genomic offset (ssp26)"), 
                                   offset_plot(data, variable = offset85, 
                                               cgr_rev = FALSE, abs_cols = T, brks = 6,
                                               plot_title = "Genomic offset (ssp85)")),
                                   ncol = 2, align = "vh")

ggsave("../plots/offset_26_85.tiff", width = 14, height = 8, dpi = 300)


# Bioclim ----------------------------------------------------------------------

bioclim <- read.csv("../data/ch_bioclim.csv") %>% 
  mutate(data = case_when(is.na(ssp) ~ "Present",
                          ssp == 85  ~ "ssp85",
                          ssp == 26  ~ "ssp26")) %>% 
  select(c("Site", "bio1", "bio5", "bio15", "data")) %>% 
  filter(!Site %in% c("Harrison", "Nahatlatch", "Raft"))

future  <- bioclim[bioclim$data == "ssp85",]   %>% select(starts_with("bio"))
present <- bioclim[bioclim$data == "Present",] %>% select(starts_with("bio"))
delta   <- cbind(bioclim[bioclim$data == "ssp85", "Site"], 
                 data.frame(c(future-present))) %>% 
                 rename("population" = 1) %>% 
                 mutate(population = gsub(" ", "", population)) %>% 
  merge(., loc, by.x = "population", by.y = "Site") 

offset <- data[,c(1:3)] %>% merge(., delta, by = "population") %>% 
  pivot_longer(cols = c("offset85", "offset26"), names_to = "offset", values_to = "offset_vals") %>% 
  pivot_longer(cols = starts_with("bio"), names_to = "biovar", values_to = "bioclim")

clab <- c(bio1 = "Annual mean temperature (bio1)",
          bio5 = "Max temperature of the warmest month (bio5)",
          bio15 = "Precipitation seasonality (bio15)",
          offset26 = "Genomic offset (ssp26)",
          offset85 = "Genomic offset (ssp85)")

(off_env_lms <- ggplot(data = offset,
       aes(x = bioclim, y = offset_vals)) +
  geom_smooth(aes(x = bioclim, y = offset_vals),
              method = "lm", colour = "black") +  theme_bw() +
  geom_point(shape = 21, aes(fill = as.factor(Latitude)),
             size = 2.5) + 
  facet_grid(offset ~ biovar, 
             scales = "free", labeller = as_labeller(clab)) +
  scale_fill_manual(values = c(viridis_pal(option = "D")(length(unique(offset$population))))) +
  labs(y = "Offset value",
       x = "Environmental change (Future - Present)") +
  stat_poly_eq(use_label(c("R2", "p")),
               label.x = "left",
               label.y = "top") +
  scale_y_continuous(expand = expansion(mult = c(1/50, 1/10))) +
  theme(legend.position = "none"))

ggsave("../plots/offset_env_var.tiff", dpi = 300,
       width = 10, height = 7)



