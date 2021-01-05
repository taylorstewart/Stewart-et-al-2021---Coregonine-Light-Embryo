#### CLEAR THE ENVIRONMENT FIRST -----------------------------------------------------------------

rm(list = ls(all.names = TRUE))


#### LOAD PACKAGES -------------------------------------------------------------------------------

library(dplyr)
library(ggplot2)


defi.files <- list.files('data/DEFI-L', full.names = TRUE, recursive = TRUE, pattern = "20191204")

## 0AA0025 - Top Shelf
## 0AA0038 - Middle Shelf
## 0AQH006 - Middle Shelf
## 0AA0037 - Bottom Shelf
## 0AA0026 - Bottom Shelf

data <- do.call(rbind, lapply(defi.files, function(i) {
  
  read.csv(i, header = TRUE, skip = 25)  %>% 
    mutate(shelf = substr(i, 30, 32),
           rep =  substr(i, 33, 33),
           timestamp = as.POSIXct(TimeStamp, format = "%Y/%m/%d %H:%M:%S")) %>% 
    select(timestamp, shelf, rep, quantum = "Quantum..umol..m.2s..")
}))

trt <- data.frame(shelf = c("Bot", "Mid", "Top"),
                  light = c("High", "Medium", "Low"))

data.time <- data %>% left_join(trt) %>% 
  mutate(time = strftime(.$timestamp, format = "%H:%M:%S"), 
         date = strftime(.$timestamp, format = "%Y-%m-%d")) %>% 
  filter(date >= "2019-12-10", time >= "11:00:00", time <= "13:00:00")


data.summary <- data.time %>% group_by(shelf, light) %>% 
  summarize(mean.quantum = mean(quantum),
            sd.quantum = sd(quantum))

