library(sp); library(terra)
library(geodata); library(tidyverse)
library(ggplot2); library(ggrepel)
library(raster)

# Set wd.
setwd("~/ots_landscape_genetics/bioclim")

# Read in site information.
sites <- read.delim("ch2023_sequenced.txt")[,c(4,8:9L)]


# ------------------------------------------------------------------------------


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
                           res  = 1/2,
                           path = "./wcdata",
                           download = FALSE)

  # Download projected bioclimatic variables.
  # Following Tigano et al. 2023 [Evo. App.]
  # Climate data at RCP2.6 (~best case).
  # Set download = TRUE for first time.
  f.clim26 <- cmip6_tile(ssp  = "126",
                         lon  = lon,
                         lat  = lat,
                         res  = 5,
                         path = "./fclim26",
                         time = "2041-2060",
                         var  = "bioc",
                         model = "UKESM1-0-LL",
                         download = FALSE)
  
  # Similar to above, but RCP8.5 (~worst case).
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
  
  # Set up a function to extract data from each of the 
  # three climate datasets isolated above.
  # Also renames columns consistently and adds two qualifiers.
  f <- function(clim, ssp, time) cbind(sites[i, 1],
       terra::extract(clim, y = points)[i, 2:20]) %>% 
       `colnames<-`(., c("Site", paste0("bio", seq(1, 19, 1)))) %>% 
       mutate(period = as.factor(time), ssp = as.factor(ssp)) 
  
  # For each site, bind together and specify time period and ssp.
  # NOTE: First (1970s-2000) data set is in a different order!!!
        # bio1, bio10, bio11,...bio9. Unlike the other two. 
        # Arrange separately otherwise there is a severe mismatch.
  c.list[[i]] <- rbind(
    f(f.clim26, time = "2041-2060", ssp = "26"),
    f(f.clim85, time = "2041-2060", ssp = "85"),
    cbind(sites[i, 1], terra::extract(c.clim, points)[i, 2:20]) %>% 
      `colnames<-`(., c("Site", paste0("bio", c(1, seq(10, 19, 1),
                        seq(2, 9, 1))))) %>% 
      mutate(period = "1970-2000", ssp = NA),
    make.row.names = FALSE # Just use numbers.
  )
}

# Put the bioclim data into a single dataframe. 
biovar <- bind_rows(c.list); head(biovar)

# No missing data.
sum(is.na(biovar)) == 0

write.csv(biovar, "../data/ch_bioclim.csv", row.names = FALSE)

# Visualize potential before/after differences.
ggplot(data = biovar %>% 
         pivot_longer(cols = paste0("bio", seq(1, 19, 1))),
       aes(x = reorder(period, value), 
           y = value, fill = ssp)) +
  scale_fill_manual(
    values = c("blue", "red"),
    limits = c("26", "85")) +
  geom_boxplot(alpha = 1/2) +
  labs(x = NULL, y = "Value") +
  facet_wrap(~name, scales = "free_y") +
  theme_bw() +
  theme(legend.position = "top",
        panel.grid = element_blank()) +
  guides(fill = guide_legend(title = "Shared Socio-economic Pathway"))


ggsave("plots/clim_ts.tiff", dpi = 300, width = 10, height = 6)


# Scaling and autocorrelation --------------------------------------------------


"%ni%" <- Negate("%in%") # Custom not-in operator.
# Retain only contemporary data and remove unnecessary columns.
cbio <- biovar[biovar$period == "1970-2000", colnames(biovar) %ni% c("period", "ssp")]

# Correlation matrix of bioclimatic variables. 
(bcor <- abs(cor(cbio[,-1]))) # Excludes site info.
bcor[!lower.tri(bcor)] <- 0; hist(bcor)

# Retain variables with less than 80% correlation.
(ev <- cbio[, !apply(bcor, 2, function(x) any(x > 0.80))])
dim(ev) # number of variables retained

# Scale and center (Z-transform) uncorrelated variables.
sval <- as.data.frame(scale(ev)) %>%  # Make numeric too.
  mutate_if(is.character, is.numeric)

# Make a population-specific dataframe with scores along various PC axes. 
(evpca <- prcomp(sval))
(pc1ve <- summary(evpca)$importance[2,1]) # PC1 % variation explained.
(pc2ve <- summary(evpca)$importance[2,2]) # PC2 % variation explained.

# Retain PC scores and attach to each site.
pop_climPC <- as.data.frame(prcomp(sval)$x) %>% 
  mutate(pop = tools::toTitleCase(tolower(gsub("\\_.*", "", cbio$Site))))

# Plot PC1 and 2 scores.
ggplot(data = pop_climPC,
       aes(x = PC1, y = PC2, label = pop)) +
  geom_point(size = 2) +
  geom_label_repel(max.overlaps = Inf,
                   min.segment.length = 0,
                   box.padding = 1/3,
                   size = 1.7,
                   segment.size = 1/5) +
  labs(y = paste0("PC2 (", round(pc2ve*100, 1), "%)"),
       x = paste0("PC1 (", round(pc1ve*100, 1), "%)")) +
  theme_bw()

ggsave("plots/envPCA.tiff", dpi = 300, width = 7, height = 7)


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
   dplyr::select(c(4,2)))

write.csv(elev, "../data/pop_elevations.csv", row.names = F)

ggplot(data = elev, 
       aes(x = reorder(pop, elevation),
           y = elevation)) + 
  geom_bar(stat = "identity",
           color = "black") +
  labs(x = NULL, y = "Elevation (m)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 40, hjust = 0.9))

ggsave("plots/elevations.tiff", dpi = 300, width = 10, height = 5)

