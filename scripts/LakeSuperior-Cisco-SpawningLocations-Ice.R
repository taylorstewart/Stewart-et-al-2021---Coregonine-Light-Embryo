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

ice.daily <- read.csv("data/SpawningLocations_Ice/LakeSuperior-DailyIce.csv", header = TRUE)


#### LOAD LAKE SUPERIOR SHAPEFILE ----------------------------------------------------------------

ls_poly <- readOGR(dsn = "data/SpawningLocations_Ice/shapefiles/lake-superior.shp", layer = "lake-superior")
ls_poly <- spTransform(ls_poly, CRS("+proj=longlat"))
ls_poly.fort <- fortify(ls_poly) %>% mutate(hole = case_when(group == 0.2 ~ "ZZZ",
                                                             group == 0.3 ~ "ZZZ",
                                                             hole == TRUE ~ "Land",
                                                             hole == FALSE ~ "Lake"))

ls_crs <- st_crs(ls_poly)


#### LOAD CISCO SPAWNING LOCATIONS ---------------------------------------------------------------

spawn <- read.csv("data/SpawningLocations_Ice/Coregonines_Superior_Spawning_Locations.csv", header = TRUE) %>% 
  filter(SPECIES == "Cisco") %>% 
  mutate(location = 1:n()) %>% 
  dplyr::select(location, lat = LATITUDE, lon = LONGITUDE)

spawn.points <- st_as_sf(spawn, coords = c("lon", "lat"), crs = ls_crs)


#### CALCULATE DATE OF 15% ICE COVER -------------------------------------------------------------

ice.files <- ice.daily %>%
  group_by(ice.year) %>% 
  filter(jday >= 46, jday < 100) %>% 
  filter(row_number() == 1) %>%
  mutate(date = gsub("-", "", date)) %>% pull(date)


#### CALCULATE NEAREST ICE NEIGHBOR FOR EACH SPAWNING LOCATION -----------------------------------

spawn.ice.all <- do.call(rbind, mclapply(ice.files, mc.cores = 8, function(i) {

  ## Load ice data for each year
  ice.data <- fread(file = paste0("/Volumes/home/GL-Seasonal-Environmental-Change/data/lake-superior-ice-concentration/daily-point/lake-superior-ice-concentration-", i, ".csv")) %>% 
    mutate(jday = yday(date),
           year = year(date),
           ice.conc = ifelse(ice.conc < 0, 0, ice.conc),
           ice.conc = ifelse(ice.conc > 100, 100, ice.conc),
           ice.year = ifelse(jday < 300, year-1, year)) %>% 
    dplyr::select(year, date, jday, ice.year, ice.conc, lon, lat)
  
  ## Create spatial points
  ice.points <- st_as_sf(ice.data, coords = c("lon", "lat"), crs = ls_crs)

  ## Find nearest ice point to each spawning location
  closest.ice <- st_nn(spawn.points, ice.points)
  closest.ice <- do.call(rbind, closest.ice)
  
  ## Extract each nearest ice value for all spawning locations
  spawn.ice <- ice.data[c(closest.ice),] %>% bind_cols(spawn) %>% 
    dplyr::select(location, ice.year, year, date, jday, ice.conc, lon = lon...6, lat = lat...7)
}))

## Create a dataframe with spawning ground coordinates
spawn.coord <- spawn.ice.all %>% filter(ice.year == 2019) %>% 
  select(location, lon, lat)

## Calculate percent
spawn.ice.perc <- spawn.ice.all %>% 
  mutate(ice.group = ifelse(ice.year < 2000, 1980, 2000)) %>% 
  group_by(ice.group) %>% 
  mutate(ice.logical = ifelse(ice.conc >= 15, 1, 0)) %>% 
  group_by(location, ice.group) %>% 
  summarize(annual.ice.perc = mean(ice.logical)) %>% 
  mutate(annual.ice.perc = ifelse(annual.ice.perc == 0, 0.1, annual.ice.perc),
         perc.bin = cut(annual.ice.perc, breaks = c(0, 0.25, 0.5, 0.75, 1),
                        labels = c("< 25%", "25-50%", "50-75%", "> 75%")),
         ice.group = factor(ice.group, ordered = TRUE,
                            levels = c(1980, 2000), labels = c("1980 - 2000", "2000 - 2020"))) %>% 
  left_join(spawn.coord)

spawn.ice.annual <- spawn.ice.all %>% select(-lat, -lon) %>% 
  left_join(spawn.coord)

# VISUALIZATION -------------------------------------------------------------------------------

ggplot(data = ls_poly.fort, aes(long, lat)) +
  geom_polygon(aes(group = group, fill = hole), 
               color = "black", size = 0.5, show.legend = FALSE) +
  scale_fill_manual(values = c("#e9f3f8", "white", "white")) + 
  new_scale_fill() + 
  geom_point(data = spawn.ice.perc, aes(x = lon, y = lat, fill = perc.bin), 
             color = "black", shape = 21, size = 3.75) +
  scale_fill_manual(values = c("#a1dab4", "#41b6c4", "#2c7fb8", "#253494")) +
  coord_fixed(ratio = 1.4) +
  scale_x_continuous(limits = c(-92.3, -84.3), breaks = seq(-92, -84, 1), expand = c(0, 0),
                     labels =  paste0(seq(-92, -84, 1), "°")) +
  scale_y_continuous(limits = c(46.35, 49.1), breaks = seq(46.5, 49, 0.5), expand = c(0, 0),
                     labels =  paste0(seq(46.5, 49, 0.5), "°")) +
  labs(x = "Longitude", y = "Latitude") +
  guides(fill = guide_legend(nrow = 1, override.aes = list(size = 6))) +
  theme(axis.text = element_text(size = 13),
        axis.title = element_blank(),
        axis.ticks.length = unit(1.5, 'mm'),
        legend.title = element_blank(), 
        legend.text = element_text(size = 13),
        legend.key = element_rect(fill = "white"),
        legend.position = c(0.233, 0.955),
        legend.background = element_rect(fill = "white", color = "black", size = 0.5, linetype = "solid"),
        strip.text = element_text(size = 14),
        strip.background = element_rect(color = "white", fill = "white"),
        panel.background = element_rect(fill = "transparent"), 
        panel.grid = element_line(linetype = 'dashed', color = "#B0B0B080"), 
        panel.border = element_rect(linetype = 'solid', color = "black", fill = "transparent"),
        panel.ontop = TRUE) +
  facet_wrap(~ice.group, nrow = 2)

ggsave("figures/LakeSuperior-Historical-Ice-CiscoSpawning.png", dpi = 300, width = 10, height = 10)


lapply(unique(spawn.ice.annual$ice.year), function(k) {
  spawn.ice.annual.filt <- spawn.ice.annual %>% filter(ice.year == k)
  
  ggplot(data = ls_poly.fort, aes(long, lat)) +
    geom_polygon(aes(group = group, fill = hole), 
                 color = "black", size = 0.5, show.legend = FALSE) +
    scale_fill_manual(values = c("#e9f3f8", "white", "white")) + 
    new_scale_fill() + 
    geom_point(data = spawn.ice.annual.filt, aes(x = lon, y = lat, fill = ice.conc), 
               color = "black", shape = 21, size = 3.75) +
    scale_fill_gradient(low = "blue", high = "white", limits = c(0.0, 100.0), breaks = c(0.0, 25.0, 50.0, 75.0, 100.0)) +
    coord_fixed(ratio = 1.4) +
    scale_x_continuous(limits = c(-92.3, -84.3), breaks = seq(-92, -84, 1), expand = c(0, 0),
                       labels =  paste0(seq(-92, -84, 1), "°")) +
    scale_y_continuous(limits = c(46.35, 49.1), breaks = seq(46.5, 49, 0.5), expand = c(0, 0),
                       labels =  paste0(seq(46.5, 49, 0.5), "°")) +
    labs(x = "Longitude", y = "Latitude") +
    guides(fill = guide_colourbar(title = paste0("15 Feb. ", k, " Ice Concentration (%)"), title.position = "top", direction = "horizontal",
                                  barheight = 1.15, barwidth = 16, 
                                  ticks.colour = "black", ticks.linewidth = 1,
                                  frame.colour = 'black', frame.linewidth = 1)) +
    theme(axis.text = element_text(size = 12),
          axis.title = element_blank(),
          axis.ticks.length = unit(1.5, 'mm'),
          legend.title = element_text(size = 13), 
          legend.text = element_text(size = 13),
          legend.key = element_rect(fill = "white"),
          legend.position = c(0.2, 0.88),
          legend.background = element_rect(fill = "white", color = "black", size = 0.5, linetype = "solid"),
          legend.margin = margin(2, 5.7, 1.5, 3, unit = 'mm'),
          panel.background = element_rect(fill = "transparent"), 
          panel.grid = element_line(linetype = 'dashed', color = "#B0B0B080"), 
          panel.border = element_rect(linetype = 'solid', color = "black", fill = "transparent"),
          panel.ontop = TRUE)
  
  ggsave(paste0("figures/Annual_Ice_Spawning/2-15-Feb/LakeSuperior-Historical-Ice-CiscoSpawning-", k, ".png"), dpi = 300, width = 10, height = 5)
})

