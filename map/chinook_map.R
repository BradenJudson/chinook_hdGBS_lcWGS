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
library(ggrepel); library(rnaturalearth)

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
sites <- read.delim(file = "ch2023_sequenced.txt") %>% 
  arrange(Latitude) %>% 
  mutate(site = tools::toTitleCase(tolower(gsub("\\_.*", "", Population))),
         sitenum = rownames(.))

# Read in global river database and subset for relevant systems.
rivers <- ne_load(scale = 10, type = "rivers_lake_centerlines", destdir = "maps", returnclass = "sf") %>% 
  filter(name %ni% c("Kuskokwim", "North Fork Kuskokwim", "S. Fork Kuskokwim", "Copper", "Kobuk", "Thelon",
                     "Susitna", "Iliamna Lake Outlet", "Coppermine", "Humboldt", "Pit", "Sacramento",
                     "Arctic Red", "Peel", "Bow", "Oldman", "Athabasca", "North Saskatchewan", "Brazeau"))

# Again, use bcmaps to download watercourse dat and subset.
(riv <- bcmaps::watercourses_5M())
(sRiv <- riv %>% filter(name_en %in% c("Peace River", "Teslin River", "Cowichan River",
                                       "Fraser River", "Lillooet River", "Adams River",
                                       "San Juan River", "Cariboo River", "Iskut River",
                                       "Taku River", "Bella Coola River", "Klinaklini River",
                                       "Taseko River", "Dean River", "Kitimat River",
                                       "Kitlope River", "Tatshenshini River", "Shuswap River",
                                       "Nass River", "Similkameen River", "Skeena River",
                                       "North Thompson River", "Thompson River", "Chilko River",
                                       "Inklin River", "Quesnel River", "Chilcotin River",
                                       "Harrison River", "Portland Canal", "Nechako River",
                                       "Stuart River", "South Thompson River")))

# Above command misses the Canadian Okanagan, so I manually read that in.
okanagan <- st_read("ok_path.kml")      

# Read in lake data downloaded from iMap BC.
# Filter some lakes. Object ID part is awkward as some lakes aren't named.
lakes <- st_transform(st_read(dsn = "DBM_BC_7H_MIL_DRAINAGE_POLY"), crs = 4326) %>% 
  filter(NAME %in% c("Harrison Lake", "Okanagan Lake", "Lower Arrow Lake", "Fran?ois Lake",
  "Upper Arrow Lake", "Kinbasket Lake", "Adams Lake", "Shuswap Lake", "Fraser Lake") | 
   OBJECTID %in% c("259", "280", "263", "260", "275", "258", "267", "225"))

(pnw <- ggplot() +
  geom_sf(data = USA, fill = "gray90", linewidth = 1/10) +
  geom_sf(data = bcn, fill = "gray90", linewidth = 1/10) +
  geom_sf(data = bch, fill = "gray90", linewidth = 1/10) +
  ggspatial::annotation_north_arrow(location = "bl",
                                    pad_x = unit(-1/10, "cm"),
             style = ggspatial::north_arrow_fancy_orienteering()) +
  ggspatial::annotation_scale(location = "bl",
                              pad_y = unit(1/10, "cm"),
                              width_hint = 1/10) + 
    geom_sf(data = okanagan, color  = "skyblue", linewidth = 1/4) +
    geom_sf(data = sRiv,     colour = "skyblue", linewidth = 1/4) +
    geom_sf(data = rivers,   colour = "skyblue", linewidth = 1/4) +  
    geom_sf(data = lakes,    colour = "skyblue", fill = "skyblue") +
    # coord_sf(xlim = c(-116, -166), ylim = c(41, 66)) +
    geom_point(data = sites, size = 2.5, stroke = 1/3,
               shape = 21, color = "black", fill = "white",
               aes(x = Longitude, y = Latitude)) +
    geom_text(data = sites, size = 1.3, fontface = "bold",
              aes(x = Longitude, y = Latitude, label = `sitenum`)) +
    scale_fill_manual(values = rep("white", nrow(sites)),
                      labels = paste(sites$sitenum, sites$site)) +
    guides(fill = guide_legend()) +
    coord_sf(xlim = c(-115, -165), ylim = c(41, 66)) +
    theme_minimal() +
    theme(legend.position = "top", panel.grid = element_blank(), legend.key = element_blank(),
          panel.background = element_rect(fill = alpha("skyblue", 1/10)),
          panel.border = element_rect(color = "black", fill = NA)))

ggsave("plots/bc_map.tiff", dpi = 300, width = 6, height = 6)

# https://stackoverflow.com/questions/24801987/numbered-point-labels-plus-a-legend-in-a-scatterplot


# Overlay CU info --------------------------------------------------------------

# From: https://open.canada.ca/data/en/dataset/2f4bd945-f47e-47e3-9108-79f6ee39242c
# CUs <- st_transform(st_read(dsn = "Chinook_Salmon_CU_Boundary", 
#                     layer = "CK_CU_BOUNDARY_En")[,1], crs = 4326); plot(CUs)

(pnwCUs <- pnw +
    # geom_sf(data = CUs, 
    #         aes(fill = CU_NAME, alpha = 1/10),
    #         linewidth = 1/20, color = "black") +
    geom_point(data = sites, size = 1/2,
               aes(x = Longitude, y = Latitude)) +
    # Right-adjusted labels.
    geom_text_repel(data = sites[sites$site %in% c("Duteau", "Adams", 
                           "Seymour", "Blue", "Slim", "Swift", 
                           "Raft", "Cariboo", "Taseko", "Okanagan", 
                           "Spius", "Nahatlatch"),], 
                    size = 1.5,
                    min.segment.length = 0, 
                    max.overlaps = Inf, 
                    nudge_y = 2,
                    segment.size = 1/10,    
                    segment.alpha = 1/2, 
                    box.padding = 0.1,
                    hjust = 1, 
                    direction = "y", 
                    xlim = c(-117.5, NA),
                    aes(label = site, 
                        x = Longitude, 
                        y = Latitude)) +
    # Left-adjusted labels.
    geom_text_repel(data = sites[sites$site %in% c("Cowichan", "Phillips",
                           "Qualicum", "Nimpkish", "Sarita", "Kaouk",
                           "Moyeha", "San", "Klinaklini", "Serpentine"),], 
                    size = 1.5, 
                    box.padding = 0.15, 
                    min.segment.length = 0, 
                    max.overlaps = Inf, 
                    segment.size = 1/10,    
                    segment.alpha = 1/2,
                    hjust = 0, 
                    direction = "y", 
                    xlim = c(NA, -129),
                    nudge_y = -3,
                    aes(label = site, 
                        x = Longitude, 
                        y = Latitude)) +
    # South-central labels.
    geom_text_repel(data = sites[sites$site %in% c("Pitt", "Cheakamus", 
                           "Birkenhead", "Portage"),], 
                    size = 1.5, 
                    nudge_y = -2, 
                    box.padding = 0.25,
                    min.segment.length = 0, 
                    max.overlaps = Inf, 
                    segment.size = 1/10,    
                    segment.alpha = 1/2,
                    direction = "both", 
                    ylim = c(49, 45.2), 
                    xlim = c(-123, NA),
                    aes(label = site, 
                        x = Longitude, 
                        y = Latitude)) +
    # North/central coast left-aligned.
    geom_text_repel(data = sites[sites$site %in% c("Yakoun", "Dean", "Nusatsum", "Kitwanga",
                           "Kildala", "Ecstall", "Kitlope", "Kilbella", "Kitsumkalum", 
                           "Kincolith", "Kwinageese", "Tahltan", "Verrett", "Takhanne"),], 
                    size = 1.5, 
                    nudge_y = -0.10, 
                    box.padding = 0.1,
                    min.segment.length = 0, 
                    max.overlaps = Inf, 
                    segment.size = 1/10,    
                    segment.alpha = 1/2,
                    direction = "y", 
                    xlim = c(NA, -140),
                    ylim = c(NA, 59),
                    aes(label = site, 
                        x = Longitude, 
                        y = Latitude)) +
    # Northern pops. 
    geom_text_repel(data = sites[sites$site %in% c(
                                                   "Big", "Hoole", "Mayo",
                                                   "Nordenskiold", "Tincup", "Klondike",
                                                    "Takhini", "Salmon", 
                                                   "Tozitna", "Andreafsky", "Trinity", 
                                                   "Imnaha", "Abernathy", "Siuslaw"),],
                    size = 1.5, nudge_x = 0.2, min.segment.length = 0, 
                    max.overlaps = Inf, nudge_y = -0.1,
                    segment.size = 1/10,    
                    segment.alpha = 1/2,
                    aes(label = site, 
                        x = Longitude, 
                        y = Latitude)) +
    # Northern pops, aligned in BC center.
    geom_text_repel(data = sites[sites$site %in% c("Sustut", "Dudidontu", "Yeth", 
                                                   "Morley", "Endako"),],
                    size = 1.5, nudge_x = 0.2, min.segment.length = 0, 
                    hjust = 1, box.padding = 0.1,
                    direction = "y", 
                    xlim = c(-124.5, NA),
                    ylim = c(56, 60),
                    segment.size = 1/10,    
                    segment.alpha = 1/2,
                    aes(label = site, 
                        x = Longitude, 
                        y = Latitude)) +
    coord_sf(xlim = c(-115, -165), ylim = c(41, 66)) +
    theme_minimal() +
    theme(legend.position = "none", panel.grid = element_blank(),
          panel.background = element_rect(fill = alpha("skyblue", 1/10)),
          panel.border = element_rect(color = "black", fill = NA)))


ggsave("plots/CUmap.tiff", dpi = 300, width = 5, height = 5)






