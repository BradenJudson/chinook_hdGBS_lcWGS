################################################################################
#
# Author: Braden J. Judson
# Email: Braden.Judson@dfo-mpo.gc.ca
# Date created: December 08, 2023
# Date last modified: April 3, 2024
#
#
################################################################################

setwd("~/ots_landscape_genetics/map")

# devtools::install_github("rspatial/geodata")

library(ggplot2); library(tidyverse); library(raster); library(sf)
library(bcmaps); library(ggspatial); library(sp); library(geodata)
library(ggrepel); library(rnaturalearth); library(cowplot)

# Custom "not in" operator.
"%ni%" <- Negate("%in%")

# Establish base map -----------------------------------------------------------

# Download American shape data and a subset of the Canadian data.
USA <- sf::st_as_sf(geodata::gadm(country = "USA", level = 0, path = "."))
prov <- sf::st_as_sf(geodata::gadm(country = "CAN", level = 1, path = "."))
bcn <- prov[prov$NAME_1 %in% c("Alberta", "Yukon", "Northwest Territories", "Nunavut"),]

# Use bcmaps package to get high resolution BC shape data.
# Necessary due to complexity of coastline.
bch <- st_transform(bc_bound_hres(), crs = 4326)

# Read in site information and reformat labels.
sites <- read.delim(file = "../data/ch2023_sequenced.txt") %>% 
  arrange(Latitude) %>% 
  mutate(site = tools::toTitleCase(tolower(gsub("\\_.*", "", Population))),
         sitenum = as.numeric(rownames(.)))

sites[sites$Population == "SAN_JUAN_RIVER", "site"] <- "San Juan"

# To make computations faster, only use data from focal area.
bound <- c(ymin = 35, ymax = 70, xmin = -170, xmax = -105)

# Data from: https://www.hydrosheds.org/products/hydrorivers
# North American and "Arctic" data are separate. 
arctic  <- st_crop(st_read("./hydroRIVERS_data/arctic_v10/HydroRIVERS_v10_ar.shp"),        y = bound)
northam <- st_crop(st_read("./hydroRIVERS_data/north_america_v10/HydroRIVERS_v10_na.shp"), y = bound)
rivers <- rbind(arctic, northam) # Combine data sources here.
rm(arctic); rm(northam); gc()    # Only retain combined object. Clears up memory.

rivers3 <- st_transform(sf::st_simplify(rivers[rivers$ORD_STRA >  2,]), crs = 4326)

# For labeling later - sets text values along y axis at predefined heights.
(yv <- 55 - (55 - 40.25)*(1:(nrow(sites)/2))/nrow(sites)*2)


(pnw <- ggplot(data = sites) +
    geom_sf(data = USA, fill = "gray90", linewidth = 1/10) +
    geom_sf(data = bcn, fill = "gray90", linewidth = 1/10) +
    geom_sf(data = bch, fill = "gray90", linewidth = 1/10) +
    ggspatial::annotation_north_arrow(location = "bl",
                                      pad_x = unit(3.65, "cm"),
                                      style = ggspatial::north_arrow_fancy_orienteering()) +
    ggspatial::annotation_scale(location = "bl", pad_x = unit(3.75, "cm"),
                                pad_y = unit(1/10, "cm"), width_hint = 1/10) +
    geom_sf(data = rivers3,   colour = "skyblue",  linewidth = 1/4) +  
    geom_point(data = sites, size = 2.5, stroke = 1/3,
               shape = 21, color = "black", fill = "white",
               aes(x = Longitude, y = Latitude)) +
    geom_text(data = sites, size = 1.3, fontface = "bold",
              aes(x = Longitude, y = Latitude, label = `sitenum`)) +
    scale_fill_manual(values = rep("white", nrow(sites)),
                      labels = paste(sites$sitenum, sites$site)) +
    guides(fill = guide_legend(override.aes = list(alpha = 0))) +
    coord_sf(xlim = c(-115, -165), ylim = c(41, 66)) +
    theme_minimal() +
    theme(legend.position = "right", panel.grid = element_blank(), 
          panel.background = element_rect(fill = alpha("skyblue", 1/10)),
          panel.border = element_rect(color = "black", fill = NA),
          plot.margin = unit(c(0,0,0,0), "cm")) +
    geom_rect(aes(ymin = 40.01, ymax = 54.8, xmax = -153, xmin = -167.25),
              colour = "black", fill = "white", alpha = 1/6, linewidth = 1/20) +
    geom_text(data = sites[sites$sitenum <= 30, ], aes(label = paste(sitenum, ". ", site), 
                                                       x = -167, y = yv, hjust = 0), size = 2, inherit.aes = F) +
    geom_text(data = sites[sites$sitenum > 30, ], aes(label = paste(sitenum, ". ", site), 
                                                      x = -160, y = yv, hjust = 0), size = 2, inherit.aes = F))


# Inset ------------------------------------------------------------------------

# Download low-res country outlines.
ca <- map_data("world", "Canada")
us <- map_data("world", "USA") 
me <- map_data("world", "Mexico")

# Make the inset plot by itself. 
(ins <- ggplot() +
    geom_polygon(data = us, aes(x = long, y = lat, group = group),
                 fill = "grey80", colour = "black", linewidth = 1/8) +
    geom_polygon(data = ca, aes(x = long, y = lat, group = group),
                 fill = "grey90", colour = "black", linewidth = 1/8) +
    geom_polygon(data = me, aes(x = long, y = lat, group = group),
                 fill = "grey70", colour = "black", linewidth = 1/8) +
    theme_void() +
    theme(panel.border = element_rect(colour = "black", 
                         fill = NA, linewidth = 1/4),
          panel.background = element_rect(fill = "white"))  +
    annotate("rect", fill = NA, colour = "black",
             linewidth = 1/2,
             xmin = -165, xmax = -115,
             ymin = 41, ymax = 66) +
    # Important to maintain accurate proportions/orientations. 
    # Plot is cartesian otherwise and appears distorted.
    coord_map(ylim = c(72, 20),
              xlim = c(-57, -165)))


# Add inset --------------------------------------------------------------------


ggdraw(plot = pnw) +
  draw_plot({
    ins
  },
  x = 0.790,
  y = 0.575,
  width = 0.2,
  height = 0.5)

ggsave("../plots/map_winset.tiff", dpi = 300, 
       width = 6, height = 6)

