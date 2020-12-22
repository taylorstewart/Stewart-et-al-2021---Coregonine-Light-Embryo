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
  dplyr::select(location, depth, lon = lon...6, lat = lat...5)


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
  dplyr::select(location, depth, lon = lon...6, lat = lat...5)


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
    dplyr::select(location, depth, ice.year, year, date, jday, ice.conc, lon = lon...6, lat = lat...7)
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
    dplyr::select(location, depth, ice.year, year, date, jday, ice.conc, lon = lon...6, lat = lat...7)
}))

ggplot(ls.spawn.ice.perc, aes(x = depth, y = ice.conc, group = factor(year), color = factor(year))) + 
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) + 
  theme_bw()


## Calculate percent
ls.spawn.ice.perc <- ls.spawn.ice.all %>% 
  filter(ice.year %in% c(2013, 2019))

lo.spawn.ice.perc <- lo.spawn.ice.all %>% 
  filter(ice.year %in% c(2018, 2019)) %>% 
  left_join(lo.spawn.coord)
