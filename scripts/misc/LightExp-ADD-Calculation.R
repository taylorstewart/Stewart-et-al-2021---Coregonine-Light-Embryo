## ===========================================================
## Clear the environment first
## ===========================================================
rm(list = ls(all.names=TRUE))


## ===========================================================
## Load packages
## ===========================================================
library(dplyr)
library(magrittr)
library(readxl)
library(ggplot2)


## ===========================================================
## Load data
## ===========================================================
climate.temp.lo.h <- read_excel("data/HOBO/2020-Artedi-Light-HOBO-Corrected.xlsx", sheet = "LO-H") %>% 
  mutate(date = as.POSIXct(paste0(substr(datetime, 1, 4), "-", substr(datetime, 6, 7), "-", substr(datetime, 9, 10))))

climate.temp.lo.m <- read_excel("data/HOBO/2020-Artedi-Light-HOBO-Corrected.xlsx", sheet = "LO-M") %>% 
  mutate(date = as.POSIXct(paste0(substr(datetime, 1, 4), "-", substr(datetime, 6, 7), "-", substr(datetime, 9, 10))))

climate.temp.lo.l <- read_excel("data/HOBO/2020-Artedi-Light-HOBO-Corrected.xlsx", sheet = "LO-L") %>% 
  mutate(date = as.POSIXct(paste0(substr(datetime, 1, 4), "-", substr(datetime, 6, 7), "-", substr(datetime, 9, 10))))

climate.temp.ls.h <- read_excel("data/HOBO/2020-Artedi-Light-HOBO-Corrected.xlsx", sheet = "LS-H") %>% 
  mutate(date = as.POSIXct(paste0(substr(datetime, 1, 4), "-", substr(datetime, 6, 7), "-", substr(datetime, 9, 10))))

climate.temp.ls.m <- read_excel("data/HOBO/2020-Artedi-Light-HOBO-Corrected.xlsx", sheet = "LS-M") %>% 
  mutate(date = as.POSIXct(paste0(substr(datetime, 1, 4), "-", substr(datetime, 6, 7), "-", substr(datetime, 9, 10))))

climate.temp.ls.l <- read_excel("data/HOBO/2020-Artedi-Light-HOBO-Corrected.xlsx", sheet = "LS-L") %>% 
  mutate(date = as.POSIXct(paste0(substr(datetime, 1, 4), "-", substr(datetime, 6, 7), "-", substr(datetime, 9, 10))))

## Combine all treatments into one dataframe
climate.temp.all <- bind_rows(climate.temp.lo.h, climate.temp.lo.m, climate.temp.lo.l,
                              climate.temp.ls.h, climate.temp.ls.m, climate.temp.ls.l)


## ===========================================================
## Calculate accumulative degree days
## ===========================================================
climate.temp.summary <- climate.temp.all %>% group_by(population, treatment, date) %>% 
  summarize(mean.daily.temp = mean(temp_c)) %>% ungroup() %>% 
  group_by(population, treatment) %>% 
  mutate(ADD = cumsum(mean.daily.temp)) %>% ungroup()

write.csv(climate.temp.summary, "data/2020-Artedi-Light-ADD.csv", row.names = FALSE)


climate.temp.mean <- climate.temp.all %>% group_by(population, treatment) %>% 
  summarize(mean.temp = mean(temp_c))


