################################################################################
#
# Author: Braden J. Judson
# Email: Braden.Judson@dfo-mpo.gc.ca
# Date created: December 08, 2023
# Date last modified: January 29, 2024
# Script description: 
#
#
################################################################################

setwd("~/ots_landscape_genetics/map")

# devtools::install_github("rspatial/geodata")

library(ggplot2); library(tidyverse); library(raster); library(sf)
library(bcmaps); library(ggspatial); library(sp); library(geodata)

"%ni%" <- Negate("%in%")


# Establish base map -----------------------------------------------------------

USA <- sf::st_as_sf(geodata::gadm(country = "USA", level = 0, path = "."))
prov <- sf::st_as_sf(geodata::gadm(country = "CAN", level = 1, path = "."))
bcn <- prov[prov$NAME_1 %in% c("Alberta", "Yukon", "Northwest Territories"),]

bch <- st_transform(bc_bound_hres(), crs = 4326)

(pnw <- ggplot() +
  geom_sf(data = USA, fill = "gray70", linewidth = 1/10) +
  geom_sf(data = bcn, fill = "gray70", linewidth = 1/10) +
  geom_sf(data = bch, fill = "gray99", linewidth = 1/10) +
  ggspatial::annotation_north_arrow(location = "bl",
              pad_x = unit(0.1, "cm"),
             style = ggspatial::north_arrow_fancy_orienteering()) +
  ggspatial::annotation_scale(location = "bl",pad_y = unit(0.1, "cm"),
                              width_hint = 1/10) + 
    coord_sf(xlim = c(-114, -139), ylim = c(48.5, 60)) +
    theme_void())

(riv <- bcmaps::watercourses_5M()) ### NEEDS FILTERING ###
# Add lakes somehow? # Or add 15M files whenever relevant....?

pnw +
  geom_sf(data = riv, colour = "skyblue")

ggsave("plots/bc_map.tiff", dpi = 300, width = 5, height = 5)



# Overlay CU info --------------------------------------------------------------

# From: https://open.canada.ca/data/en/dataset/2f4bd945-f47e-47e3-9108-79f6ee39242c
CUs <- st_transform(st_read(dsn = "Chinook_Salmon_CU_Boundary", 
                    layer = "CK_CU_BOUNDARY_En")[,1], crs = 4326); plot(CUs)

sites <- read.delim(file = "ch2023_sequenced.txt")

(pnwCUs <- pnw +
    geom_point(data = sites, aes(x = Longitude, y = Latitude)) +
    geom_sf(data = CUs, 
            aes(fill = CU_NAME, alpha = 1/10),
            linewidth = 1/20, color = "black") +
    coord_sf(xlim = c(-115, -165),
             ylim = c(41, 66)) +
    theme_minimal() +
    theme(legend.position = "none", panel.grid = element_blank(),
          panel.background = element_rect(fill = alpha("skyblue", 1/10)),
          panel.border = element_rect(color = "black", fill = NA)))


ggsave("plots/CUmap.tiff", dpi = 300, width = 5, height = 5)





