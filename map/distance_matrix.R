setwd("~/ots_landscape_genetics/map")

library(sf); library(ggplot2); library(tidyverse); library(rmapshaper)


# Resources ---------------------------------------------------------------

# https://ggplot2.tidyverse.org/reference/ggsf.html#combining-sf-layers-and-regular-geoms
# https://stackoverflow.com/questions/28595414/calculate-distance-between-2-lon-lats-but-avoid-going-through-a-coastline-in-r
# https://gis.stackexchange.com/questions/353633/r-spatial-erase-one-polygon-from-another-correct-use-of-st-difference
# https://osf.io/an6b5/download [use st_difference?]
# https://gis.stackexchange.com/questions/109639/reverse-clipping-erasing-in-r
# https://stackoverflow.com/questions/22947448/finding-euclidean-distance-in-rspatstat-between-points-confined-by-an-irregul
# https://stackoverflow.com/questions/55893463/measuring-the-distance-of-all-points-along-a-line-in-r-linestring-in-sf

# st_buffer returns polygon?

# Make polygon of pacific ocean basically with this:
# https://gis.stackexchange.com/questions/403977/sf-create-polygon-from-minimum-x-and-y-coordinates
# Remove overlap of land, getting ocean and rivers
# follow link above to compute distances
# adjust lat/lon as needed
# Consider river order 3 or 4? 

# -------------------------------------------------------------------------

# To make computations faster, only use data from focal area.
bound <- c(ymin = 40, ymax = 70, xmin = -170, xmax = -115)

# Data from: https://www.hydrosheds.org/products/hydrorivers
# North American and "Arctic" data are separate. 
arctic  <- st_crop(st_read("./hydroRIVERS_data/arctic_v10/HydroRIVERS_v10_ar.shp"),        y = bound)
northam <- st_crop(st_read("./hydroRIVERS_data/north_america_v10/HydroRIVERS_v10_na.shp"), y = bound)
rivers <- rbind(arctic, northam) # Combine data sources here.
rm(arctic); rm(northam); gc()    # Only retain combined object. Clears up memory.
  

# Transform each linear feature into a polygon by expanding it with a constant buffer. 
# Benefits from some geometric simplification (faster & alternative is unnecessarily complex).
# Other transformations are necessary for st_difference later on.
buff <- sf::st_buffer(sf::st_simplify(rivers[rivers$ORD_STRA > 3,]), dist = 1500)
water_poly  <- st_make_valid(st_union(buff)) # Collapse overlapping geometries.
water_multi <- st_crop(st_union(st_make_valid(st_cast(water_poly,   # Reformat.
                       to = "MULTIPOLYGON"))), y = bound)           # Crop.
st_write(water_multi, "./hydroRIVERS_data/water_union.shp", append = FALSE) # Write in case things go wrong.


# Polygons for Canada and USA. Combine and remove subsets. 
USA <- sf::st_as_sf(geodata::gadm(country = "USA", level = 0, path = "."))
CAN <- sf::st_as_sf(geodata::gadm(country = "CAN", level = 0, path = "."))
NAM <- st_crop(rbind(USA, CAN), y = bound); plot(NAM[,1]); rm(USA); rm(CAN)


# Compute a polygon that is just the land, not the ocean or rivers.
land <- st_difference(x = NAM, y = water_multi, dimension = "polygon")
ggplot(NULL) + geom_sf(data = land, fill = "white", colour = NA) +
  theme_classic() + theme(panel.background = element_rect(fill = 'gray90')) 
ggsave("../plots/land_only.png", width  = 1000, height = 1000, units = "px")


rr <- raster(extent(land), nrow = 50, ncol = 50)
j <- rasterize(land, rr)

j[is.na(j)] <- 0

tr1 <- geoCorrection(transition(j, transitionFunction = mean, directions = 16), type = "c")








