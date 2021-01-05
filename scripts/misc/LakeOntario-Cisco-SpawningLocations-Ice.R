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


#### LOAD ICE DATA -------------------------------------------------------------------------------

ice.daily <- read.csv("data/SpawningLocations_Ice/LakeOntario-DailyIce.csv", header = TRUE)


#### LOAD LAKE SUPERIOR SHAPEFILE ----------------------------------------------------------------

lo.poly <- readOGR(dsn = "data/SpawningLocations_Ice/shapefiles/lake-ontario.shp", layer = "lake-ontario")
lo.poly <- spTransform(lo.poly, CRS("+proj=longlat"))
lo.poly.fort <- fortify(lo.poly)

lo.crs <- st_crs(lo.poly)


#### LOAD CISCO SPAWNING LOCATIONS ---------------------------------------------------------------

lo.spawn <- read.csv("data/SpawningLocations_Ice/Coregonines_Ontario_Spawning_Locations.csv", header = TRUE) %>% 
  dplyr::mutate(location = 1:n()) %>% 
  dplyr::select(location, lat, lon)

lo.spawn.points <- st_as_sf(lo.spawn, coords = c("lon", "lat"), crs = lo.crs)


#### CALCULATE DATE OF 15% ICE COVER -------------------------------------------------------------

ice.files <- ice.daily %>%
  group_by(ice.year) %>% 
  filter(jday >= 15, jday < 100) %>% 
  filter(row_number() == 1) %>%
  mutate(date = gsub("-", "", date)) %>% pull(date)


#### CALCULATE NEAREST ICE NEIGHBOR FOR EACH SPAWNING LOCATION -----------------------------------

lo.spawn.ice.all <- do.call(rbind, mclapply(ice.files, mc.cores = 8, function(i) {
  
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
lo.spawn.coord <- lo.spawn.ice.all %>% filter(ice.year == 2019) %>% 
  select(location, lon, lat)

## Calculate percent
lo.spawn.ice.perc <- lo.spawn.ice.all %>% 
  filter(ice.year %in% c(2017, 2018)) %>% 
  #mutate(ice.group = case_when(ice.year >= 1980 & ice.year < 1990 ~ "1980 - 1990",
  #                             ice.year >= 2010 & ice.year <= 2020 ~ "2010 - 2020",
  #                             TRUE ~ "remove")) %>% 
  #mutate(ice.group = case_when(ice.year >= 1980 & ice.year < 1990 ~ "1980 - 1990",
  #                             ice.year >= 1990 & ice.year < 2000 ~ "1990 - 2000",
  #                             ice.year >= 2000 & ice.year < 2010 ~ "2000 - 2010",
  #                             ice.year >= 2010 & ice.year <= 2020 ~ "2010 - 2020",
  #                             TRUE ~ "remove")) %>% 
  #filter(ice.group != "remove") %>% 
  #group_by(ice.group, location) %>% 
  group_by(year, location) %>% 
  mutate(ice.logical = ifelse(ice.conc >= 15, 1, 0)) %>% 
  #dplyr::summarize(annual.ice.perc = mean(ice.logical)) %>% 
  #mutate(annual.ice.perc = ifelse(annual.ice.perc == 0, 0.1, annual.ice.perc),
  #       perc.bin = cut(annual.ice.perc, breaks = c(0, 0.25, 0.5, 0.75, 1),
  #                      labels = c("< 25%", "25-50%", "50-75%", "> 75%"))) %>% #,
  #ice.group = factor(ice.group, ordered = TRUE)) %>% 
  left_join(lo.spawn.coord)

lo.spawn.ice.annual <- lo.spawn.ice.all %>% select(-lat, -lon) %>% 
  left_join(lo.spawn.coord)


#### VISUALIZATION -------------------------------------------------------------------------------

ggplot(data = lo.poly.fort, aes(long, lat)) +
  geom_polygon(aes(group = group, fill = hole), 
               color = "black", size = 0.5, show.legend = FALSE) +
  scale_fill_manual(values = c("#e9f3f8", "white", "white")) + 
  new_scale_fill() + 
  geom_point(data = lo.spawn.ice.perc, aes(x = lon, y = lat, fill = ice.conc), 
             color = "black", shape = 21, size = 3.75) +
  scale_fill_gradient(low = "blue", high = "white", limits = c(0.0, 100.0), breaks = c(0.0, 25.0, 50.0, 75.0, 100.0)) +
  coord_fixed(ratio = 1.4) +
  scale_x_continuous(limits = c(-80.0, -75.85), breaks = seq(-80, -75, 0.5), expand = c(0, 0),
                     labels =  paste0(seq(-80, -75, 0.5), "°")) +
  scale_y_continuous(limits = c(43.12, 44.45), breaks = seq(43.25, 44.75, 0.25), expand = c(0, 0),
                     labels =  paste0(seq(43.25, 44.75, 0.25), "°")) +
  labs(x = "Longitude", y = "Latitude") +
  guides(fill = guide_colourbar(title = "Ice Concentration (%)", title.position = "top", direction = "horizontal",
                                barheight = 1.15, barwidth = 16, 
                                ticks.colour = "black", ticks.linewidth = 1,
                                frame.colour = 'black', frame.linewidth = 1)) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_blank(),
        axis.ticks.length = unit(1.5, 'mm'),
        legend.title = element_text(size = 13), 
        legend.text = element_text(size = 13),
        legend.key = element_rect(fill = "white"),
        legend.position = c(0.2, 0.94),  ## Fill Bar Legend
        legend.background = element_rect(fill = "white", color = "black", size = 0.5, linetype = "solid"),
        legend.margin = margin(2, 5.7, 1.5, 3, unit = 'mm'),
        strip.text = element_text(size = 14),
        strip.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"), 
        panel.grid = element_line(linetype = 'dashed', color = "#B0B0B080"), 
        panel.border = element_rect(linetype = 'solid', color = "black", fill = "transparent"),
        panel.ontop = TRUE) +
  facet_wrap(~year, nrow = 2)

ggsave("figures/LakeOntario-Historical-Ice-CiscoSpawning-FillBar.png", dpi = 300, width = 10, height = 9.7)

