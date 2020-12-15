#### CLEAR THE ENVIRONMENT FIRST -----------------------------------------------------------------

rm(list = ls(all.names = TRUE))


#### LOAD PACKAGES -------------------------------------------------------------------------------

library(readxl)
library(dplyr)
library(magrittr)
library(ggplot2)
library(lubridate)
library(tidyr)


#### LOAD LIGHT DATA -----------------------------------------------------------------------------

light <- read_excel('data/Lake_Superior_Ice_Light/DEFI2-L_20161026_LakeSuperior.xlsx', sheet = "DEFI2-L") %>% 
  mutate(timestamp = ymd_hms(TimeStamp)) %>%   ## Sensors deployed in Eastern TZ
  select(-TimeStamp)

## Filter by dates deployed
light %<>% filter(timestamp >= ymd_hms("2016-11-30 01:00:00"), timestamp <= ymd_hms("2017-06-02 01:00:00"))

## Combine years
light.filt <- light %>% 
  mutate(timestamp = as.POSIXct(timestamp, format = "%Y-%m-%d %H:%M:%S"),
         date = as.POSIXct(strptime(timestamp, format = "%Y-%m-%d")),
         time = as.POSIXct(strftime(timestamp, format = "%H:%M:%S", tz = "UTC"), format = "%H:%M:%S"),
         year = year(date),
         jday = yday(date),
         ice.year = ifelse(jday > 300, year+1, year)) %>% 
  filter(!is.na(time)) %>% 
  select(ice.year, date, time, light = Quantum, -timestamp, -Batt)

## Filter to daylight hours
light.filt.day <- light.filt %>% filter(time >= as.POSIXct(paste0(Sys.Date(), " 09:00:00")), time <= as.POSIXct(paste0(Sys.Date(), " 15:00:00")))


#### LOAD ICE DATA -------------------------------------------------------------------------------

## Create list of file names (by year)
ice.files <- list.files('data/Lake_Superior_Ice_Light/Sand-Island-Ice', full.names = TRUE, pattern = ".xls")

## Read in all files and combine into one dataset.
ice.full <- do.call(bind_rows, lapply(ice.files, function(i) {
  print(i)
  read_excel(i, sheet = paste0("T", substr(i, 46, 53))) %>% 
    mutate(date = as.POSIXct(paste0(substr(i, 46, 49), "-", substr(i, 50, 51), "-", substr(i, 52, 53))),
           year = year(date),
           jday = yday(date),
           ice.year = ifelse(jday > 300, year+1, year),
           CT = as.numeric(ifelse(is.na(CT) == TRUE, 00, CT))) %>% 
    dplyr::select(date, ice.conc = CT)
}))


#### COMBINE LIGHT AND ICE DATA IN TIDYR FORMAT --------------------------------------------------

ice.light.tidy <- left_join(light.filt, ice.full) %>% 
  mutate(ice.conc = ifelse(is.na(ice.conc) == TRUE, 0, ice.conc)) %>% 
  pivot_longer(4:5, names_to = "measurement", values_to = "value") %>% 
  mutate(timestamp = as.POSIXct(paste0(date, " ", format(time, format = "%H:%M:%S"))))


#### COMBINE DAY LIGHT AND ICE DATA --------------------------------------------------------------

ice.light.day <- left_join(light.filt.day, ice.full) %>% 
  mutate(ice.conc = ifelse(is.na(ice.conc) == TRUE, 0, ice.conc))


#### FILTER TO EACH ICE TREATMENT ----------------------------------------------------------------

icelight.10 <- ice.light.day %>% filter(ice.conc < 10) %>% 
  group_by(ice.year) %>% 
  summarize(n = n(),
            mean.light = mean(light),
            sd.light = sd(light),
            mean.ice = mean(ice.conc))

icelight.40.60 <- ice.light.day %>% filter(ice.conc > 40, ice.conc < 60) %>% 
  group_by(ice.year) %>% 
  summarize(n = n(),
            mean.light = mean(light),
            sd.light = sd(light),
            mean.ice = mean(ice.conc))

icelight.90 <- ice.light.day %>% filter(ice.conc > 90) %>% 
  group_by(ice.year) %>% 
  summarize(n = n(),
            mean.light = mean(light),
            sd.light = sd(light),
            mean.ice = mean(ice.conc))


## ===========================================================
## Plot time-series
## ===========================================================
ggplot(ice.light.tidy, aes(x = timestamp, y = value, group = measurement)) +
  geom_path(aes(colour = measurement), size = 1.5) +
  scale_x_datetime(date_labels = "%m/%Y", date_breaks = "1 month", expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(10, 100, 10), expand = c(0, 0)) +
  scale_color_manual(values = c("#0091e6", "gray50"), labels = c("Ice", "Light")) +
  labs(y = expression(paste("Ice Coverage (%)\n\nPhoton Flux (μmol ", m^-2, " ", s^-1, ")", sep="")), x = "") +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25), 
        axis.ticks.length = unit(2, 'mm'),
        legend.position = c(0.08, 0.92), 
        legend.text = element_text(size = 20),
        legend.title = element_blank(), 
        legend.key.width = unit(2.5, 'lines'), 
        legend.key.height = unit(1.5, 'lines'),
        legend.key = element_rect(fill = "white"), 
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        plot.margin = unit(c(5, 20, 0, 20), "mm"))

ggsave("figures/2017-Ice-Light.png", dpi = 300, width = 15, height = 7.5)
