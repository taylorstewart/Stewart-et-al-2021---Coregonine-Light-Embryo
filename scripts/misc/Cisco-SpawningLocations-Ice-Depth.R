#### CLEAR THE ENVIRONMENT FIRST -----------------------------------------------------------------

rm(list = ls(all.names = TRUE))


#### LOAD PACKAGES -------------------------------------------------------------------------------

library(dplyr)
library(data.table)
library(ggplot2)
library(rgdal)
library(ggpolypath)
library(sf)
library(nngeo)
library(ggnewscale)
library(parallel)
library(gridExtra)
library(grid)
library(cowplot)
library(raster)


#### LOAD ICE DATA -------------------------------------------------------------------------------

ls.ice.daily <- read.csv("data/SpawningLocations_Ice/LakeSuperior-DailyIce.csv", header = TRUE)
lo.ice.daily <- read.csv("data/SpawningLocations_Ice/LakeOntario-DailyIce.csv", header = TRUE)


#### LOAD LAKE SUPERIOR SHAPEFILE ----------------------------------------------------------------

ls.poly <- readOGR(dsn = "data/SpawningLocations_Ice/shapefiles/lake-superior.shp", layer = "lake-superior")
ls.poly <- spTransform(ls.poly, CRS("+proj=longlat"))
ls.poly.fort <- fortify(ls.poly) %>% mutate(hole = case_when(group == 0.2 ~ "ZZZ",
                                                             group == 0.3 ~ "ZZZ",
                                                             hole == TRUE ~ "Land",
                                                             hole == FALSE ~ "Lake"))
ls.crs <- st_crs(ls.poly)


#### LOAD LAKE SUPERIOR SHAPEFILE ----------------------------------------------------------------

lo.poly <- readOGR(dsn = "data/SpawningLocations_Ice/shapefiles/lake-ontario.shp", layer = "lake-ontario")
lo.poly <- spTransform(lo.poly, CRS("+proj=longlat"))
lo.poly.fort <- fortify(lo.poly)
lo.crs <- st_crs(lo.poly)


#### LOAD CISCO SPAWNING LOCATIONS ---------------------------------------------------------------

ls.spawn <- read.csv("data/SpawningLocations_Ice/Coregonines_Superior_Spawning_Locations.csv", header = TRUE) %>% 
  filter(SPECIES == "Cisco") %>% 
  dplyr::mutate(location = 1:n()) %>% 
  dplyr::select(location, lat = LATITUDE, lon = LONGITUDE)
lo.spawn <- read.csv("data/SpawningLocations_Ice/Coregonines_Ontario_Spawning_Locations.csv", header = TRUE) %>% 
  dplyr::mutate(location = 1:n()) %>% 
  dplyr::select(location, lat, lon)

ls.spawn.points <- st_as_sf(ls.spawn, coords = c("lon", "lat"), crs = ls.crs)
lo.spawn.points <- st_as_sf(lo.spawn, coords = c("lon", "lat"), crs = lo.crs)


#### LAKE SUPERIOR BATHYMETRY RASTER TO POINTS ---------------------------------------------------

ls.bathy <- readOGR(dsn = "data/SpawningLocations_Ice/shapefiles/lake-superior-bathymetry-point.shp", layer = "lake-superior-bathymetry-point")
ls.bathy.sf <- st_as_sf(ls.bathy, coords = c("lon", "lat"), crs = ls.crs)
ls.bathy.df <- as(ls.bathy, 'data.frame') %>% 
  dplyr::select(depth = grid_code, lon = coords.x1, lat = coords.x2)

## Find nearest ice point to each spawning location
closest.ice <- st_nn(ls.spawn.points, ls.bathy.sf)
closest.ice <- do.call(rbind, closest.ice)

## Extract depth and bind to spawning locations
ls.spawn.depth <- ls.bathy.df[c(closest.ice),] %>% dplyr::bind_cols(ls.spawn) %>% 
  dplyr::select(location, depth, lon = lon...6, lat = lat...5) %>% 
  mutate(depth = as.numeric(gsub("-", "", depth)),
         depth.cat = case_when(depth < 20 ~ "0-20",
                               depth >= 20 & depth < 40 ~ "20-40",
                               depth >= 40 & depth < 60 ~ "40-60",
                               depth >= 60 & depth < 80 ~ "60-80",
                               depth >= 80 & depth <= 100 ~ ">80"),
         depth.cat = factor(depth.cat, levels = c("0-20", "20-40", "40-60", "60-80", ">80")))
  

# LAKE ONTARIO BATHYMETRY RASTER TO POINTS ----------------------------------------------------

lo.bathy <- readOGR(dsn = "data/SpawningLocations_Ice/shapefiles/lake-ontario-bathymetry-point.shp", layer = "lake-ontario-bathymetry-point")
lo.bathy <- spTransform(lo.bathy, CRS("+proj=longlat"))
lo.bathy.sf <- st_as_sf(lo.bathy, coords = c("lon", "lat"), crs = lo.crs)
lo.bathy.df <- as(lo.bathy, 'data.frame') %>% 
  dplyr::select(depth = grid_code, lon = coords.x1, lat = coords.x2)

## Find nearest ice point to each spawning location
lo.closest.ice <- st_nn(lo.spawn.points, lo.bathy.sf)
lo.closest.ice <- do.call(rbind, lo.closest.ice)

## Extract depth and bind to spawning locations
lo.spawn.depth <- lo.bathy.df[c(lo.closest.ice),] %>% dplyr::bind_cols(lo.spawn) %>% 
  dplyr::select(location, depth, lon = lon...6, lat = lat...5) %>% 
  mutate(depth = as.numeric(gsub("-", "", depth)),
         depth.cat = case_when(depth < 5 ~ "0-5",
                               depth >= 5 & depth < 10 ~ "5-10",
                               depth >= 10 & depth < 15 ~ "10-15",
                               depth >= 15 & depth < 20 ~ "15-20",
                               depth >= 20 & depth < 25 ~ "20-25"))


#### CALCULATE DATE OF 15% ICE COVER -------------------------------------------------------------

ls.ice.files <- ls.ice.daily %>%
  group_by(ice.year) %>% 
  filter(jday >= 31, jday < 100) %>% 
  filter(row_number() == 1) %>%
  mutate(date = gsub("-", "", date)) %>% pull(date)

lo.ice.files <- lo.ice.daily %>%
  group_by(ice.year) %>% 
  filter(jday >= 31, jday < 100) %>% 
  filter(row_number() == 1) %>%
  mutate(date = gsub("-", "", date)) %>% pull(date)


#### CALCULATE NEAREST ICE NEIGHBOR FOR EACH SPAWNING LOCATION -----------------------------------

ls.spawn.ice.all <- do.call(rbind, mclapply(ls.ice.files, mc.cores = 8, function(i) {
  
  ## Load ice data for each year
  ice.data <- fread(file = paste0("/Volumes/home/GL-Seasonal-Environmental-Change/data/lake-superior-ice-concentration/daily-point/lake-superior-ice-concentration-", i, ".csv")) %>% 
    dplyr::mutate(jday = yday(date),
                  year = year(date),
                  ice.conc = ifelse(ice.conc < 0, 0, ice.conc),
                  ice.conc = ifelse(ice.conc > 100, 100, ice.conc),
                  ice.year = ifelse(jday < 300, year-1, year)) %>% 
    dplyr::select(year, date, jday, ice.year, ice.conc, lon, lat)
  
  ## Create spatial points
  ice.points <- st_as_sf(ice.data, coords = c("lon", "lat"), crs = ls.crs)
  
  ## Find nearest ice point to each spawning location
  closest.ice <- st_nn(ls.spawn.points, ice.points)
  closest.ice <- do.call(rbind, closest.ice)
  
  ## Extract each nearest ice value for all spawning locations
  spawn.ice <- ice.data[c(closest.ice),] %>% dplyr::bind_cols(ls.spawn.depth) %>% 
    dplyr::select(location, depth, depth.cat, ice.year, year, date, jday, ice.conc, lon = lon...6, lat = lat...7)
}))

lo.spawn.ice.all <- do.call(rbind, mclapply(lo.ice.files, mc.cores = 8, function(i) {
  
  ## Load ice data for each year
  ice.data <- fread(file = paste0("/Volumes/home/GL-Seasonal-Environmental-Change/data/lake-ontario-ice-concentration/lake-ontario-ice-concentration-", i, ".csv")) %>% 
    dplyr::mutate(jday = yday(date),
                  year = year(date),
                  ice.conc = ifelse(ice.conc < 0, 0, ice.conc),
                  ice.conc = ifelse(ice.conc > 100, 100, ice.conc),
                  ice.year = ifelse(jday < 300, year-1, year)) %>% 
    dplyr::select(year, date, jday, ice.year, ice.conc, lon, lat)
  
  ## Create spatial points
  ice.points <- st_as_sf(ice.data, coords = c("lon", "lat"), crs = lo.crs)
  
  ## Find nearest ice point to each spawning location
  closest.ice <- st_nn(lo.spawn.points, ice.points)
  closest.ice <- do.call(rbind, closest.ice)
  
  ## Extract each nearest ice value for all spawning locations
  spawn.ice <- ice.data[c(closest.ice),] %>% dplyr::bind_cols(lo.spawn.depth) %>% 
    dplyr::select(location, depth,  depth.cat, ice.year, year, date, jday, ice.conc, lon = lon...6, lat = lat...7)
}))

## Create a dataframe with spawning ground coordinates
ls.spawn.coord <- ls.spawn.ice.all %>% filter(ice.year == 2019) %>% 
  dplyr::select(location, depth, depth.cat, lon, lat)
lo.spawn.coord <- lo.spawn.ice.all %>% filter(ice.year == 2019) %>% 
  dplyr::select(location, depth, depth.cat, lon, lat)

## Calculate percent
ls.spawn.ice.coord <- ls.spawn.ice.all %>% 
  left_join(ls.spawn.coord) %>% 
  filter(year >= 2000)

lo.spawn.ice.coord <- lo.spawn.ice.all %>% 
  left_join(lo.spawn.coord) %>% 
  filter(year >= 2000)


#### VISUALIZATIONS ------------------------------------------------------------------------------

## Lake Superior
plot.ls.map <- ggplot(data = ls.poly.fort, aes(long, lat)) +
  geom_polygon(aes(group = group, fill = hole), 
               color = "gray50", size = 0.5, show.legend = FALSE) +
  scale_fill_manual(values = c("#e9f3f8", "white", "white")) + 
  new_scale_fill() + 
  geom_point(data = ls.spawn.depth, aes(x = lon, y = lat, fill = depth.cat), 
             color = "black", shape = 21, size = 3.75) +
  scale_fill_manual(values = c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c")) +
  coord_fixed(ratio = 1.4) +
  scale_x_continuous(limits = c(-92.3, -84.3), breaks = seq(-92, -84, 1), expand = c(0, 0),
                     labels =  paste0(seq(-92, -84, 1), "°")) +
  scale_y_continuous(limits = c(46.35, 49.1), breaks = seq(46.5, 49, 0.5), expand = c(0, 0),
                     labels =  paste0(seq(46.5, 49, 0.5), "°")) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_blank(),
        axis.ticks.length = unit(1.5, 'mm'),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.position = "top",
        legend.title = element_blank(),
        legend.key = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"), 
        panel.grid = element_line(linetype = 'dashed', color = "#B0B0B080"), 
        panel.border = element_rect(linetype = 'solid', color = "black", fill = "transparent"),
        panel.ontop = TRUE)

plot.ls.heatmap <- ggplot(ls.spawn.ice.coord, aes(x = year, y = location, group = depth.cat, fill = ice.conc)) +
  geom_tile() + 
  scale_fill_gradient(low = "blue", high = "white", limits = c(0.0, 100.0), breaks = c(0.0, 25.0, 50.0, 75.0, 100.0)) +
  scale_x_continuous(breaks = seq(2000, 2020, 1), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "", y = "Spawning Locations") +
  guides(fill = guide_colourbar(title = "Ice Concentration (%)", title.position = "top", direction = "horizontal",
                                barheight = 1, barwidth = 15, 
                                ticks.colour = "black", ticks.linewidth = 1,
                                frame.colour = 'black', frame.linewidth = 1)) +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_blank(),
        axis.ticks.length.y = unit(0, "mm"),
        legend.text = element_text(size = 13),
        legend.position = "top") +
  facet_grid(depth.cat ~ ., space = "free", scales = "free")

## Lake Ontario
plot.lo.map <- ggplot(data = lo.poly.fort, aes(long, lat)) +
  geom_polygon(aes(group = group, fill = hole), 
               color = "gray50", size = 0.5, show.legend = FALSE) +
  scale_fill_manual(values = c("#e9f3f8", "white", "white")) + 
  new_scale_fill() + 
  geom_point(data = lo.spawn.depth, aes(x = lon, y = lat, fill = depth.cat), 
             color = "black", shape = 21, size = 3.75) +
  scale_fill_manual(values = c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c")) +
  coord_fixed(ratio = 1.5) +
  scale_x_continuous(limits = c(-80.0, -75.85), breaks = seq(-80, -75, 0.5), expand = c(0, 0),
                     labels =  paste0(seq(-80, -75, 0.5), "°")) +
  scale_y_continuous(limits = c(43.12, 44.45), breaks = seq(43.25, 44.75, 0.25), expand = c(0, 0),
                     labels =  paste0(seq(43.25, 44.75, 0.25), "°")) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_blank(),
        axis.ticks.length = unit(1.5, 'mm'),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.position = "top",
        legend.title = element_blank(),
        legend.key = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"), 
        panel.grid = element_line(linetype = 'dashed', color = "#B0B0B080"), 
        panel.border = element_rect(linetype = 'solid', color = "black", fill = "transparent"),
        panel.ontop = TRUE)

plot.lo.heatmap <- ggplot(lo.spawn.ice.coord, aes(x = year, y = location, group = depth.cat, fill = ice.conc)) +
  geom_tile() + 
  scale_fill_gradient(low = "blue", high = "white", limits = c(0.0, 100.0), breaks = c(0.0, 25.0, 50.0, 75.0, 100.0)) +
  scale_x_continuous(breaks = seq(2000, 2020, 1), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "", y = "Spawning Locations") +
  guides(fill = guide_colourbar(title = "Ice Concentration (%)", title.position = "top", direction = "horizontal",
                                barheight = 1, barwidth = 15, 
                                ticks.colour = "black", ticks.linewidth = 1,
                                frame.colour = 'black', frame.linewidth = 1)) +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_blank(),
        axis.ticks.length.y = unit(0, "mm"),
        legend.text = element_text(size = 13),
        legend.position = "top") +
  facet_grid(depth.cat ~ ., space = "free", scales = "free")

## Combine all
plot.all <- grid.arrange(
  arrangeGrob(
    plot.ls.map,
    plot.lo.map,
    plot.ls.heatmap,
    plot.lo.heatmap,
    ncol = 2,
    nrow = 2,
    heights = c(1, 2.5),
    widths = c(1, 1)
  )
)

ggsave("figures/Historical-Ice-CiscoSpawning-FillBar.png", plot = plot.all, dpi = 300, width = 20, height = 18)


