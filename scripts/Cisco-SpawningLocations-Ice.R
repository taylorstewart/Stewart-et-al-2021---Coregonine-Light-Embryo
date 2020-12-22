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
  spawn.ice <- ice.data[c(closest.ice),] %>% dplyr::bind_cols(ls.spawn) %>% 
    dplyr::select(location, ice.year, year, date, jday, ice.conc, lon = lon...6, lat = lat...7)
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
    dplyr::select(location, ice.year, year, date, jday, ice.conc, lon = lon...6, lat = lat...7)
}))

## Create a dataframe with spawning ground coordinates
ls.spawn.coord <- ls.spawn.ice.all %>% filter(ice.year == 2019) %>% 
  select(location, lon, lat)
lo.spawn.coord <- lo.spawn.ice.all %>% filter(ice.year == 2019) %>% 
  select(location, lon, lat)

## Calculate percent
ls.spawn.ice.perc <- ls.spawn.ice.all %>% 
  filter(ice.year %in% c(2018, 2019)) %>% 
  left_join(ls.spawn.coord)

ls.mean.ice <- ls.spawn.ice.perc %>% 
  group_by(year) %>% 
  summarize(mean.ice = paste0(round(mean(ice.conc), 1), " %"))


lo.spawn.ice.perc <- lo.spawn.ice.all %>% 
  filter(ice.year %in% c(2018, 2019)) %>% 
  left_join(lo.spawn.coord)

lo.mean.ice <- lo.spawn.ice.perc %>% 
  group_by(year) %>% 
  summarize(mean.ice = paste0(round(mean(ice.conc), 1), " %"))


#### VISUALIZATION -------------------------------------------------------------------------------

plot.ls.2019 <- ggplot(data = ls.poly.fort, aes(long, lat)) +
  geom_polygon(aes(group = group, fill = hole), 
               color = "gray50", size = 0.5, show.legend = FALSE) +
  scale_fill_manual(values = c("#e9f3f8", "white", "white")) + 
  new_scale_fill() + 
  geom_point(data = filter(ls.spawn.ice.perc, year == 2019), aes(x = lon, y = lat, fill = ice.conc), 
             color = "black", shape = 21, size = 3.75) +
  scale_fill_gradient(low = "blue", high = "white", limits = c(0.0, 100.0), breaks = c(0.0, 25.0, 50.0, 75.0, 100.0)) +
  coord_fixed(ratio = 1.4) +
  scale_x_continuous(limits = c(-92.3, -84.3), breaks = seq(-92, -84, 1), expand = c(0, 0),
                     labels =  paste0(seq(-92, -84, 1), "°")) +
  scale_y_continuous(limits = c(46.35, 49.1), breaks = seq(46.5, 49, 0.5), expand = c(0, 0),
                     labels =  paste0(seq(46.5, 49, 0.5), "°")) +
  annotate("text", label = paste0("Mean Ice Concentration: ", filter(ls.mean.ice, year == 2019)$mean.ice), 
           x = -92.15, y = 48.375, hjust = 0, size = 5) + 
  labs(x = "Longitude", y = "Latitude", title = "31 January 2019") +
  guides(fill = guide_colourbar(title = "Ice Concentration (%)", title.position = "top", direction = "horizontal",
                                barheight = 1, barwidth = 15, 
                                ticks.colour = "black", ticks.linewidth = 1,
                                frame.colour = 'black', frame.linewidth = 1)) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_blank(),
        axis.ticks.length = unit(1.5, 'mm'),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 13), 
        legend.text = element_text(size = 13),
        legend.key = element_rect(fill = "white"),
        legend.position = c(0.197, 0.88),  ## Fill Bar Legend
        legend.background = element_rect(fill = "white", color = "black", size = 0.5, linetype = "solid"),
        legend.margin = margin(2, 5.7, 1.5, 3, unit = 'mm'),
        panel.background = element_rect(fill = "transparent"), 
        panel.grid = element_line(linetype = 'dashed', color = "#B0B0B080"), 
        panel.border = element_rect(linetype = 'solid', color = "black", fill = "transparent"),
        panel.ontop = TRUE)

plot.ls.2020 <- ggplot(data = ls.poly.fort, aes(long, lat)) +
  geom_polygon(aes(group = group, fill = hole), 
               color = "gray50", size = 0.5, show.legend = FALSE) +
  scale_fill_manual(values = c("#e9f3f8", "white", "white")) + 
  new_scale_fill() + 
  geom_point(data = filter(ls.spawn.ice.perc, year == 2020), aes(x = lon, y = lat, fill = ice.conc), 
             color = "black", shape = 21, size = 3.75) +
  scale_fill_gradient(low = "blue", high = "white", limits = c(0.0, 100.0), breaks = c(0.0, 25.0, 50.0, 75.0, 100.0)) +
  coord_fixed(ratio = 1.4) +
  scale_x_continuous(limits = c(-92.3, -84.3), breaks = seq(-92, -84, 1), expand = c(0, 0),
                     labels =  paste0(seq(-92, -84, 1), "°")) +
  scale_y_continuous(limits = c(46.35, 49.1), breaks = seq(46.5, 49, 0.5), expand = c(0, 0),
                     labels =  paste0(seq(46.5, 49, 0.5), "°")) +
  annotate("text", label = paste0("Mean Ice Concentration: ", filter(ls.mean.ice, year == 2020)$mean.ice), 
           x = -92.15, y = 48.975, hjust = 0, size = 5) + 
  labs(x = "Longitude", y = "Latitude", title = "31 January 2020") +
  guides(fill = guide_colourbar(title = "Ice Concentration (%)", title.position = "top", direction = "horizontal",
                                barheight = 1.15, barwidth = 16, 
                                ticks.colour = "black", ticks.linewidth = 1,
                                frame.colour = 'black', frame.linewidth = 1)) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_blank(),
        axis.ticks.length = unit(1.5, 'mm'),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent"), 
        panel.grid = element_line(linetype = 'dashed', color = "#B0B0B080"), 
        panel.border = element_rect(linetype = 'solid', color = "black", fill = "transparent"),
        panel.ontop = TRUE)

plot.lo.2019 <- ggplot(data = lo.poly.fort, aes(long, lat)) +
  geom_polygon(aes(group = group, fill = hole), 
               color = "gray50", size = 0.5, show.legend = FALSE) +
  scale_fill_manual(values = c("#e9f3f8", "white", "white")) + 
  new_scale_fill() + 
  geom_point(data = filter(lo.spawn.ice.perc, year == 2019), aes(x = lon, y = lat, fill = ice.conc), 
             color = "black", shape = 21, size = 3.75) +
  scale_fill_gradient(low = "blue", high = "white", limits = c(0.0, 100.0), breaks = c(0.0, 25.0, 50.0, 75.0, 100.0)) +
  coord_fixed(ratio = 1.5) +
  scale_x_continuous(limits = c(-80.0, -75.85), breaks = seq(-80, -75, 0.5), expand = c(0, 0),
                     labels =  paste0(seq(-80, -75, 0.5), "°")) +
  scale_y_continuous(limits = c(43.12, 44.45), breaks = seq(43.25, 44.75, 0.25), expand = c(0, 0),
                     labels =  paste0(seq(43.25, 44.75, 0.25), "°")) +
  annotate("text", label = paste0("Mean Ice Concentration: ", filter(lo.mean.ice, year == 2019)$mean.ice), 
           x = -79.95, y = 44.39, hjust = 0, size = 5) + 
  labs(x = "Longitude", y = "Latitude", title = "31 January 2019") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_blank(),
        axis.ticks.length = unit(1.5, 'mm'),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent"), 
        panel.grid = element_line(linetype = 'dashed', color = "#B0B0B080"), 
        panel.border = element_rect(linetype = 'solid', color = "black", fill = "transparent"),
        panel.ontop = TRUE)

plot.lo.2020 <- ggplot(data = lo.poly.fort, aes(long, lat)) +
  geom_polygon(aes(group = group, fill = hole), 
               color = "gray50", size = 0.5, show.legend = FALSE) +
  scale_fill_manual(values = c("#e9f3f8", "white", "white")) + 
  new_scale_fill() + 
  geom_point(data = filter(lo.spawn.ice.perc, year == 2020), aes(x = lon, y = lat, fill = ice.conc), 
             color = "black", shape = 21, size = 3.75) +
  scale_fill_gradient(low = "blue", high = "white", limits = c(0.0, 100.0), breaks = c(0.0, 25.0, 50.0, 75.0, 100.0)) +
  coord_fixed(ratio = 1.5) +
  scale_x_continuous(limits = c(-80.0, -75.85), breaks = seq(-80, -75, 0.5), expand = c(0, 0),
                     labels =  paste0(seq(-80, -75, 0.5), "°")) +
  scale_y_continuous(limits = c(43.12, 44.45), breaks = seq(43.25, 44.75, 0.25), expand = c(0, 0),
                     labels =  paste0(seq(43.25, 44.75, 0.25), "°")) +
  annotate("text", label = paste0("Mean Ice Concentration: ", filter(lo.mean.ice, year == 2020)$mean.ice), 
           x = -79.95, y = 44.39, hjust = 0, size = 5) + 
  labs(x = "Longitude", y = "Latitude", title = "31 January 2020") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_blank(),
        axis.ticks.length = unit(1.5, 'mm'),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent"), 
        panel.grid = element_line(linetype = 'dashed', color = "#B0B0B080"), 
        panel.border = element_rect(linetype = 'solid', color = "black", fill = "transparent"),
        panel.ontop = TRUE)


plot.all <- grid.arrange(
  arrangeGrob(
    plot.ls.2019,
    plot.lo.2019,
    plot.ls.2020,
    plot.lo.2020,
    ncol = 2,
    widths = c(1, 1)
  )
)

ggsave("figures/Historical-Ice-CiscoSpawning-FillBar.png", plot = plot.all, dpi = 300, width = 20, height = 10)

