setwd("~/ots_landscape_genetics/map")

library(sf); library(ggplot2); library(tidyverse); library(rmapshaper)
library(raster); library(gdistance); library(fasterize); library(viridis)
library(gganimate)

# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------

# To make computations faster, only use data from focal area.
bound <- c(ymin = 35, ymax = 70, xmin = -170, xmax = -115)

res <- 250000
(height <- geosphere::distm(c(bound[3], bound[1]), c(bound[3], bound[2])))
(width  <- geosphere::distm(c(bound[3], bound[1]), c(bound[4], bound[1])))
(x.pixels <- round(width/sqrt(res),  0))
(y.pixels <- round(height/sqrt(res), 0))

# Data from: https://www.hydrosheds.org/products/hydrorivers
# North American and "Arctic" data are separate. 
arctic  <- st_crop(st_read("./hydroRIVERS_data/arctic_v10/HydroRIVERS_v10_ar.shp"),        y = bound)
northam <- st_crop(st_read("./hydroRIVERS_data/north_america_v10/HydroRIVERS_v10_na.shp"), y = bound)
rivers <- rbind(arctic, northam) # Combine data sources here.
rm(arctic); rm(northam); gc()    # Only retain combined object. Clears up memory.


# Transform each linear feature into a polygon by expanding it with a constant-width buffer zone. 
# Benefits from some geometric simplification (faster & alternative is unnecessarily complex).
# Other transformations are necessary for st_difference later on.
buff <- sf::st_buffer(sf::st_simplify(rivers[rivers$ORD_STRA > 3,]), dist = 5000)
water_poly  <- st_make_valid(st_union(buff)) # Collapse overlapping geometries.
water_multi <- st_crop(st_union(st_make_valid(st_cast(water_poly,   # Reformat.
                       to = "MULTIPOLYGON"))), y = bound)           # Crop.
#st_write(water_multi, "./hydroRIVERS_data/water_union.shp", append = FALSE) # Write in case things go wrong.
water_multi <- st_read("./hydroRIVERS_data/water_union.shp")

# Polygons for Canada and USA. Combine and remove subsets. 
USA <- sf::st_as_sf(geodata::gadm(country = "USA", level = 0, path = "."))
CAN <- sf::st_as_sf(geodata::gadm(country = "CAN", level = 0, path = "."))

# Compute a polygon that is just the land, not the ocean or rivers.
# Have to do this shape-by-shape otherwise there are mysterious holes in the borders.
CANLAND <- st_difference(x = CAN, y = water_multi, dimension = "polygon")
USALAND <- st_difference(x = USA, y = water_multi, dimension = "polygon")
rm(USA); rm(CAN); gc() # Tidy up memory so I can run this locally.

# Polygon of what is essentially the Pacific Ocean in our study area.
bound <- c(ymin = 35, ymax = 70, xmin = -170, xmax = -115)
pcd <- data.frame(lon = c(-110, -170), lat = c(30, 80)) %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>% 
  st_bbox() %>% st_as_sfc() %>% st_cast("POLYGON") 

sf_use_s2(FALSE) # Required for st_difference on class MULTIPOYLGON(s).
# Remove Canadian landmass from the above polygon.
water2 <- st_difference(x = pcd, y = CANLAND, dimension = "polygon")
# Remove the American landmass from the above polygon.
water3 <- st_difference(x = water2, y = USALAND, dimension = "polygon")
# Borders become incomplete if done in a single step. 


# Simplify water geometries to reduce computational demand during rasterization.
# ms_simplify much, much faster than st_simplify. 
water_s <- ms_simplify(water3, keep = 1/50, keep_shapes = TRUE) %>% 
  st_transform(crs = '+proj=aeqd +lat_0=53.6 +lon_0=12.7') %>% 
  st_buffer(., dist = 1000) %>% st_transform(crs = 4326)
rm(water2); rm(water3); gc()

ggplot(NULL) + geom_sf(data = water_s, fill = "skyblue") +
  scale_x_continuous(limits = c(-115, -140)) +
  scale_y_continuous(limits = c(50, 60))
ggsave("../plots/simplified_waterways.png")



# -------------------------------------------------------------------------



# Make a raster of the area of interest. 
# Cols and rows are the number of cells (i.e., the resolution).
rr <- raster(extent(st_crop(USALAND, y = bound)), nrow = 5e3, ncol = 5e3*(y.pixels/x.pixels))                
rast <- fasterize(water_s, rr)  # Convert land polygon into above raster.

png("../plots/raster10000.png", width = 1000, height = 1000); plot(rast, col = '#0868ac', 
                                                              ylim = c(40, 75)); dev.off()
writeRaster(rast, "../data/invraster_10000.grd", overwrite = TRUE)
r10000 <- raster("../data/invraster_10000.grd")

sf_use_s2(TRUE) # Reinstate spherical (s2) geometry defaults.
# Prepare sparse transition matrix for all adjacent and next-to-adjacent cells (direction = 8).
# Direction 8 is only by connecting points, so it cannot skip over NA regions (i.e., jump between rivers). 
tr1 <- geoCorrection(transition(r10000, transitionFunction = mean, directions = 8), type = "c")
# saveRDS(tr1, "waterway_transition_matrix_5k.rds")
tr1 <- readRDS("waterway_transition_matrix_5k.rds")

# Read in site information.
sites <- read.delim("../data/ch2023_sequenced_adj.txt", sep = "\t")

# Convert sites to an object of class SpatialPointsDataFrame.
sitepoints <- SpatialPoints(coords = sites[,c("Longitude", "Latitude")],
                            proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

# Empty list to populate.
spLine_list <- list()

# Nested for-loop that isolates unique and uni-directional combinations of populations.
for (j in 1:nrow(sitepoints@coords)) {
  
  for (k in j:nrow(sitepoints@coords)) {
    
  if(k == j) { next } # Skip self-comparisons (i.e., distance = 0).
    
    # Display loop progres over pairs of populations.
    pops <- sites$Site
    print(paste0("Calculating distance from ", pops[j], " to ", pops[k]))
  
  # Shortest path between sites.
  SP <- shortestPath(tr1, sitepoints[j], sitepoints[k], output = "SpatialLines")
  
  # Fortify lines object for plotting purposes. 
  # # Also compute and add inter-population distances.
  SPDF <- fortify(SpatialLinesDataFrame(SP, data = data.frame(ID = 1))) %>%
    mutate(pop1 = sites[j, "Site"], pop2 = sites[k, "Site"],
           distance_m = round(geosphere::lengthLine(SP), 0))

  # Subset dataframe above and poplate list.
  spLine_list[[length(spLine_list)+1]] <- SPDF
  
  }
}

# Reformat list of spatial lines into a dataframe. 
# Make a 'label' for plotting purposes.
site_lines <- do.call("rbind", spLine_list) %>% 
  mutate(pair = as.factor(paste(pop1, "-", pop2)),
         label = as.factor(paste0(pair, ": ", format(distance_m, format = "d", big.mark = ","), " (m)")))
# write_rds(site_lines, "distance_lines.rds")
# site_lines <- readRDS("distance_lines.rds")


# Isolate unique combinations of populations to summarize distances.
dists <- site_lines %>% group_by(pair) %>%
  sample_n(1) %>% .[,7:10] # select columns 7-10.
hist(dists$distance_m, main = NULL, xlab = "Distance (m)")

# Transform inter-population distances into a pairwise matrix.
# Add column and row that are all zeroes so it can be flipped across the diagonal.
distmat <- as.data.frame.matrix(xtabs(distance_m ~ ., dists[,c(1:3)]))
distmat <- cbind(Abernathy = 0, distmat); distmat <- rbind(distmat, Yeth = 0)
distmat[4,5:39] <- distmat[5:39,4]; distmat[5:39,4] <- 0 # Fix Big Qualicum for some reason.
distmat[lower.tri(distmat)] <- t(distmat[upper.tri(distmat)]) # Make it symmetrical across the diagonal.
colnames(distmat) <- rownames(distmat) <- gsub(" ", "", rownames(distmat)) # Make naming consistent.
sum(distmat == 0) == nrow(distmat) # Check that formatting is what we expect.
write.csv(distmat, "../data/Otsh_distances_mat.csv")


# Plot distance paths between all populations. 
# Some formatting required to pass ggplot object to gganimate and have it
# function as expected (i.e., only certain layers are animated).
# Subset Takhanne as an example. Must reorder the factor after filtering.
Takhanne <- site_lines[site_lines$pop1 == "Takhanne" | site_lines$pop2 == "Takhanne",] %>% 
  arrange(distance_m) %>% mutate(label = as.factor(as.character(label)))
(pop_dists <- ggplot(data = Takhanne, 
                    aes(x = long, y = lat, colour = label, group = label)) +
    geom_sf(data = CANLAND, fill = "gray", colour = NA, inherit.aes = FALSE) +
    geom_sf(data = USALAND, fill = "gray", colour = NA, inherit.aes = FALSE) +
  # geom_raster(data = dfsp, aes(x = x, y = y, fill = layer), inherit.aes = FALSE) +
  scale_fill_manual(values = c("gray", "gray99")) +
  scale_x_continuous(expand = c(0,0), limits = c(-170, -115)) +
  scale_y_continuous(expand = c(0,0),
                     limits = c(39, 72)) + coord_sf() +
  geom_path(linewidth = 1.5) +
    scale_colour_manual(values = c(viridis_pal(option = "D")(length(unique(Takhanne$label))))) +
    geom_point(data = sites, aes(x = Longitude, y = Latitude), inherit.aes = FALSE) +
    theme_void() + theme(legend.position = "none"))
ggsave("../plots/pop_distances.tiff", width = 7, height = 8, dpi = 300)

# Define animations for ggplot object above.
# Print label in bottom left corner.
(dists_an <- pop_dists +
  transition_states(label) + 
  ease_aes('linear') +
  labs(tag = "{closest_state}") +
  theme(plot.tag.position = c(0.15, 0.1),
        plot.tag = element_text(size = 20, hjust = 0)))

# Animate based on criteria described above. 10 frames per transition and loop infinitely. 
(gif <- animate(dists_an, width = 1000, height = 1000, renderer = gifski_renderer(loop = TRUE), 
               nframes = 10*length(unique((Takhanne$pair)))))

# Remove temporary png files and save output gif.
file.remove(list.files(pattern=".png"))
anim_save("../plots/takhanne_distances.gif", gif)
