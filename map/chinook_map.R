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

setwd("~/landscape_genetics/map")

# devtools::install_github("rspatial/geodata")

library(ggplot2); library(tidyverse); library(raster); library(sf)
library(bcmaps); library(ggspatial); library(sp); library(geodata)

"%ni%" <- Negate("%in%")


# Establish base map -----------------------------------------------------------

# https://bcgov.github.io/bcmaps/articles/bcmaps.html
# https://bcgov.github.io/bcmaps/reference/index.html
# See above fo 1:5mil base map of water? 

bcn <- bc_neighbours() 
neighbours <- bcn[bcn$name %ni% c("British Columbia", "Pacific Ocean"),]
ocean <- bcn[bcn$name == "Pacific Ocean",]
bc <- bc_bound()

(pnw <- ggplot() +
  geom_sf(data = ocean, fill = "skyblue", alpha = 2/5) +
  geom_sf(data = neighbours, fill = "gray80") +
  geom_sf(data = bc, fill = "gray99") +
  ggspatial::annotation_north_arrow(location = "bl", 
                                    pad_y = unit(4/3, "cm"), 
                                    pad_x = unit(2/3, "cm"),
                                    style = ggspatial::north_arrow_fancy_orienteering()) +
  ggspatial::annotation_scale(location = "bl",
                              pad_x = unit(1, "cm"),
                              pad_y = unit(1, "cm"),
                              width_hint = 1/10) + theme_void())

(riv <- bcmaps::watercourses_5M()) ### NEEDS FILTERING ###

pnw +
  geom_sf(data = riv, colour = "skyblue")

ggsave("plots/bc_map.tiff", dpi = 300, width = 5, height = 5)



# Overlay CU info --------------------------------------------------------------

# Need to filter CUs not analyzed?
# Get Yukon/NWT shapefile. See: https://stackoverflow.com/questions/24072075/mapping-specific-states-and-provinces-in-r
# Rivers from https://www2.gov.bc.ca/gov/content/data/geographic-data-services/topographic-data/freshwater
# or using osmdata package?
# some discrepencies between coastal borders and CU boundaries.... Looks awful.
# Maybe 


CUs <- st_transform(st_read(dsn = "Chinook_Salmon_CU_Boundary", 
                    layer = "CK_CU_BOUNDARY_En")[,1], crs = 4326); plot(CUs)

(pnwCUs <- pnw +
    geom_sf(data = CUs, color = "gray", alpha = 1/2, linewidth = 1/20) +
    coord_sf(ylim = c(47, 65),
             xlim = c(-150, -113),
             crs = "WGS84"))


ggsave("plots/CUmap.tiff", dpi = 300, width = 5, height = 5)





