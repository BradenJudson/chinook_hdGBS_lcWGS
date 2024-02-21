library(sp); library(terra)
library(geodata); library(tidyverse)
library(ggplot2); library(ggrepel)

# Set wd.
setwd("~/ots_landscape_genetics/bioclim")

# Read in site information.
sites <- read.delim("ch2023_sequenced.txt")[,c(4,8:9L)]


# Retrieve bioclimatic data ----------------------------------------------------


# Establish a blank list for populating.
c.list <- list()

# For each site, extract bioclimatic data.
for (i in 1:nrow(sites)) {
  
  # Isolate site-specific lat and long.
  lon <- sites[i, "Longitude"]; lat <- sites[i, "Latitude"]
  
  # Download data from WorldClim at 5 arcsmin.
  # Store downloaded *.tiff files locally.
  # First, download current (1970 - 2000) conditions.
  c.clim <- worldclim_tile(var  = "bio",
                           lon  = lon, 
                           lat  = lat, 
                           res  = 5,
                           path = "./wcdata", 
                           download = FALSE)
  
  # Download project bioclimatic variables.
  f.clim26 <- cmip6_tile(ssp  = "126",
                         lon  = lon,
                         lat  = lat,
                         res  = 5,
                         path = "./fclim26",
                         time = "2041-2060",
                         var  = "bioc",
                         model = "UKESM1-0-LL",
                         download = FALSE)
  
  f.clim85 <- cmip6_tile(ssp  = "585",
                         lon  = lon,
                         lat  = lat,
                         res  = 5,
                         path = "./fclim85",
                         time = "2041-2060",
                         var  = "bioc",
                         model = "UKESM1-0-LL",
                         download = FALSE)

  # Extract site coordinates and project.
  points <- vect(sites[,c(2:3)], crs = "EPSG:4326",
                 geom = c("Longitude", "Latitude"))
  
  # Join with site name and populate blank list.
  c.list[[i]] <- cbind(sites[i, 1], 
  # Extract bioclim data and rename the columns.
  terra::extract(c.clim, y = points)[i,2:20]) %>% 
    `colnames<-`(., c("Site", paste0("bio", seq(1, 19, 1)))) %>% 
    mutate(grp = "current")
  
  
  
}

# Put the bioclim data into a single dataframe. 
biovar <- bind_rows(c.list); head(biovar)

# No missing data.
sum(is.na(biovar)) == 0

write.csv(biovar, "ch_bioclim.csv", row.names = FALSE)

# Scaling and autocorrelation --------------------------------------------------

# Correlation matrix of bioclimatic variables. 
(bcor <- abs(cor(biovar[,2:ncol(biovar)])))
bcor[!lower.tri(bcor)] <- 0; hist(bcor)

# Retain variables with less than 80% correlation.
(ev <- biovar[, !apply(bcor, 2, function(x) any(x > 0.80))])

# Scale and center (Z-transform) uncorrelated variables.
sval <- as.data.frame(scale(ev)) %>%  # Make numeric too.
  mutate_if(is.character, is.numeric)

# Make a population-specific dataframe with scores along various PC axes. 
pop_climPC <- as.data.frame(prcomp(sval[,1:(ncol(sval)-1)])$x) %>% 
  mutate(pop = tools::toTitleCase(tolower(gsub("\\_.*", "", biovar$Site))))

ggplot(data = pop_climPC,
       aes(x = PC1, y = PC2, label = pop)) +
  geom_point(size = 2) +
  geom_label_repel(max.overlaps = Inf,
                   min.segment.length = 0,
                   box.padding = 1/3,
                   size = 1.7,
                   segment.size = 1/5) +
  labs(y = "PC2 (9.49%)", 
       x = "PC1 (82.0%)") +
  theme_bw()

ggsave("plots/envPCA.tiff", dpi = 300, width = 5, height = 5)


# Migration harshness ----------------------------------------------------------

# Sensu Moore et al. 2017; Molecular Ecology. 

# Can measure KMLs using st_length. 
# But most/all points are at the mouth of the respective body of water.
# E.g., Cowichan appears to have a migration distance of 2m. 
# Where do they spawn? Hatchery sites? Consult w/ Tim.

# Extract site-specific elevation values and convert to dataframe.
(elev <- elevatr::get_elev_point(locations = sites[,c(3,2)] %>% 
                               `colnames<-`(., c("x", "y")),
                             units = "meters", 
                             prj = "EPSG:4326",
                             src = "aws") %>% 
   as.data.frame(.) %>% 
   mutate(pop = tools::toTitleCase(tolower(gsub("\\_.*",
                                  "", sites$Population)))) %>% 
   select(c(4,2)))

write.csv(elev, "pop_elevations.csv", row.names = F)

ggplot(data = elev, 
       aes(x = reorder(pop, elevation),
           y = elevation)) + 
  geom_bar(stat = "identity",
           color = "black") +
  labs(x = NULL, y = "Elevation (m)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 40, hjust = 0.9))

ggsave("plots/elevations.tiff", dpi = 300, width = 10, height = 5)

