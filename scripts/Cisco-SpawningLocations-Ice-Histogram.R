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
library(ggridges)


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

ls.spawn <- data.frame(lat = 46.85, lon = -90.55)
lo.spawn <- data.frame(lat = 44.05, lon = -76.20)

ls.spawn.points <- st_as_sf(ls.spawn, coords = c("lon", "lat"), crs = ls.crs)
lo.spawn.points <- st_as_sf(lo.spawn, coords = c("lon", "lat"), crs = lo.crs)


#### CALCULATE DATE OF 15% ICE COVER -------------------------------------------------------------

ls.ice.files <- ls.ice.daily %>%
  group_by(ice.year) %>% 
  filter(jday >= 1, jday < 75) %>% 
  #filter(row_number() == 1) %>%
  mutate(date = gsub("-", "", date)) %>% pull(date)

lo.ice.files <- lo.ice.daily %>%
  group_by(ice.year) %>% 
  filter(jday >= 1, jday < 75) %>% 
  #filter(row_number() == 1) %>%
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
  spawn.ice <- ice.data[c(closest.ice),] %>% dplyr::bind_cols(ls.spawn) %>% 
    dplyr::select(ice.year, year, date, jday, ice.conc, lon = lon...6, lat = lat...7)
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
  spawn.ice <- ice.data[c(closest.ice),] %>% dplyr::bind_cols(lo.spawn) %>% 
    dplyr::select(ice.year, year, date, jday, ice.conc, lon = lon...6, lat = lat...7)
}))

ls.spawn.ice.mean <- ls.spawn.ice.all %>% group_by(year) %>% 
  summarize(mean.ice.conc = mean(ice.conc))

lo.spawn.ice.mean <- lo.spawn.ice.all %>% group_by(year) %>% 
  summarize(mean.ice.conc = mean(ice.conc))

spawn.ice.all <- bind_rows(ls.spawn.ice.mean, lo.spawn.ice.mean) %>% 
  mutate(lake = factor(rep(c("Superior", "Ontario"), each = 48), ordered = TRUE,
                      levels = c("Ontario", "Superior")))


#### VISUALIZATION -------------------------------------------------------------------------------

ggplot(spawn.ice.all, aes(x = mean.ice.conc, y = lake, fill = lake, color = lake)) +
  #geom_density(alpha = 0.2) +
  #geom_jitter(aes(x = mean.ice.conc)) +
  geom_density_ridges(alpha = 0.2, point_alpha = 1,  point_shape = 16, point_size = 3,
                      jittered_points = TRUE, quantile_lines = TRUE, quantiles = 4, 
                      scale = 0.95, show.legend = FALSE) +
  scale_fill_manual(values = c("#b2df8a", "#a6cee3")) +
  scale_color_manual(values = c("#33a02c", "#1f78b4")) +
  scale_x_continuous(limit = c(0, 100), expand = c(0, 0.75)) +
  #scale_y_discrete(expand = c(0.05, 0)) +
  scale_y_discrete(expand = expansion(add = c(0.0, 0.575))) +
  labs(x = "Mean Ice Concentration (%)") +
  theme_bw() +
  theme(axis.title.x = element_text(color = "Black", size = 22, margin = margin(10, 0, 0, 0)),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 22),
        axis.ticks.length = unit(2, 'mm'),
        legend.position = "none",
        plot.margin = unit(c(5, 5, 5, 5), 'mm'))

ggsave("figures/Historical-Ice-CiscoSpawning-Histogram.png", dpi = 300, width = 10, height = 7.5)


ggplot(spawn.ice.all, aes(x = mean.ice.conc, group = lake, fill = lake)) +
  geom_histogram(bins = 30, color = "black") +
  scale_fill_manual(values = c("#b2df8a", "#a6cee3")) +
  scale_y_continuous(limits = c(0, 10), expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 1)) +
  labs(x = "Ice Concentration (%)", y = "Frequency", fill = "") +
  theme_bw() +
  theme(axis.title.x = element_text(color = "Black", size = 22, margin = margin(10, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 22, margin = margin(0, 10, 0, 0)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.ticks.length = unit(2, 'mm'),
        legend.position = "top",
        legend.text = element_text(size = 15),
        plot.margin = unit(c(5, 5, 5, 5), 'mm'))


ggplot(spawn.ice.all, aes(x = mean.ice.conc, group = lake, fill = lake)) +
  geom_density(alpha = 0.4) +
  scale_fill_manual(values = c("#b2df8a", "#a6cee3")) +
  scale_x_continuous(limits = c(0, 100), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 0.026), expand = c(0, 0)) +
  labs(x = "Mean Ice Concentration (%)", y = "Density", fill = "") +
  theme_bw() +
  theme(axis.title.x = element_text(color = "Black", size = 22, margin = margin(10, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 22, margin = margin(0, 10, 0, 0)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.ticks.length = unit(2, 'mm'),
        legend.position = "top",
        legend.text = element_text(size = 15),
        plot.margin = unit(c(5, 5, 5, 5), 'mm'))

