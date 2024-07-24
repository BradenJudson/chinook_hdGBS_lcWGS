setwd("~/ots_landscape_genetics/map")

library(sf); library(ggplot2); library(tidyverse); library(rmapshaper)
library(raster); library(gdistance)


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
bound <- c(ymin = 35, ymax = 70, xmin = -170, xmax = -115)

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


# Make a raster of the area of interest. 
# Cols and rows are the number of cells (i.e., the resolution).
rr <- raster(extent(land), nrow = 1e4, ncol = 1e4)
rast <- rasterize(land, rr)  # Convert land polygon into above raster.
invrast <- is.na(rast); plot(invrast) # Essentially invert the raster (land is NA, water is >0). 

png("../plots/raster1000.png", width = 1000, height = 1000); plot(invrast); dev.off()
writeRaster(invrast, "../data/invraster_1000.grd")
r1000 <- raster("../data/invraster_1000.grd")

# Prepare sparse transition matrix for all adjacent and next-to-adjacent cells (direction = 16).
tr1 <- geoCorrection(transition(r1000, transitionFunction = mean, directions = 16), type = "c")

sites <- read.delim("../data/ch2023_sequenced.txt") %>%
  filter(Site %in% c("Cowichan", "Harrison")) %>%
  .[,c("Site", "Longitude", "Latitude")]
coords <- sites[,c(2:3)]



# -------------------------------------------------------------------------

testdat <- data.frame(
  Site = c("a", "b", "c", "d", "e"),
  lon = c(-130, -140, -140, -123.7, -162.692),
  lat  = c(41, 50, 70, 49.3, 61.9048))

coords <- testdat[,c(2:3)]

# Convert sites to an object of class SpatialPointsDataFrame.
sitepoints <- SpatialPoints(coords = coords,
                            proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

# Empty list to populate.
spLine_list <- list()

# Nested for-loop that isolates unique and uni-directional combinations of populations.
for (j in 1:nrow(sitepoints@coords)) {
  
  for (k in j:nrow(sitepoints@coords)) {
    
  if(k == j) { next } # Skip self-comparisons (i.e., distance = 0).
  
  # Shortest path between sites.
  SP <- shortestPath(tr1, sitepoints[j], sitepoints[k], output = "SpatialLines")
  
  # Fortify lines object for plotting purposes. 
  # # Also compute and add inter-population distances.
  SPDF <- fortify(SpatialLinesDataFrame(SP, data = data.frame(ID = 1))) %>%
    mutate(pop1 = testdat[j, "Site"], pop2 = testdat[k, "Site"],
           distance_m = round(geosphere::lengthLine(SP), 0))

  # Subset dataframe above and poplate list.
  spLine_list[[length(spLine_list)+1]] <- SPDF
  
  }
}

# Reformat list of spatial lines into a dataframe. 
site_lines <- do.call("rbind", spLine_list) %>% 
  mutate(pair = as.factor(paste0(pop1, "-", pop2)))

# Plot map and shortest paths between populations. 
(dist_plot <- ggplot(NULL) + coord_sf() +
  geom_raster(data = dfsp, aes(x = x, y = y, fill = layer)) +
  geom_path(data = site_lines, aes(x = long, y = lat, colour = pair)) +
  geom_label(data = testdat, aes(x = lon, y = lat, label = Site), 
             label.r = unit(1/2, units = "lines")) +
  # geom_point(data = testdat, aes(x = lon, y = lat)) + 
  scale_fill_manual(values = c("gray", "gray99")) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0),
                     limits = c(39, 72)) +
  theme_void() + theme(legend.position = "none"))
ggsave("../plots/pop_distances.tiff", width = 7, height = 8, dpi = 300)

# Isolate unique combinations of populations to summarize distances.
dists <- site_lines %>% group_by(pair) %>%
  sample_n(1) %>% .[,7:10] # select columns 7-10.
hist(dists$distance_m, main = NULL, xlab = "Distance (m)")

# Convert distances into a pairwise matrix.
dist_mat <- dists %>% 
  pivot_wider(id_cols = pop2, 
              names_from = pop1, 
              values_from = distance_m) %>% 
  column_to_rownames("pop2")
write.csv(dist_mat, "../data/pop_water_distances.csv", row.names = T)
