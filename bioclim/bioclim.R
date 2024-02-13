library(sp); library(terra)
library(geodata); library(tidyverse)

# Set wd.
setwd("~/ots_landscape_genetics/bioclim")

# Read in site information.
sites <- read.delim("ch2023_sequenced.txt")[,c(4,8:9L)]

# Establish a blank list for populating.
c.list <- list()

# For each site, extract bioclimatic data.
for (i in 1:nrow(sites)) {
  
  # Isolate site-specific lat and long.
  lon <- sites[i, "Longitude"]; lat <- sites[i, "Latitude"]
  
  # Download data from WorldClim at 30s resolution (~1km^2).
  # Store downloaded *.tiff files locally.
  data <- worldclim_tile(var = "bio", lon = lon, 
                         lat = lat, res = 1/2,
                         path = tempdir())
  
  # Extract site coordinates and project.
  points <- vect(sites[,c(2:3)], crs = "EPSG:4326",
                 geom = c("Longitude", "Latitude"))
  
  # Join with site name and populate blank list.
  c.list[[i]] <- cbind(sites[i, 1], 
  # Extract bioclim data and rename the columns.
  terra::extract(data, y = points)[i,2:20]) %>% 
    `colnames<-`(., c("Site", paste0("bio", seq(1, 19, 1))))
  
}

# Put the bioclim data into a single dataframe. 
biovar <- bind_rows(c.list)

# Some issues with site coords being in the ocean. Ask TH.
biovar[is.na(biovar$bio1), "Site"]

