library(dplyr)
library(ggplot2)

defi.files <- list.files('/Users/Taylor/CloudStation/Cisco-Climate-Change/NA-Coregonine-Latitude/data/DEFI',
                        full.names = TRUE, recursive = TRUE, pattern = "20191204")

## 0AA0025 - Top Shelf
## 0AA0038 - Middle Shelf
## 0AQH006 - Middle Shelf
## 0AA0037 - Bottom Shelf
## 0AA0026 - Bottom Shelf

data <- do.call(rbind, lapply(defi.files, function(i) {
  
  read.csv(i, header = TRUE, skip = 25)  %>% 
    mutate(shelf = substr(i, 99, 101),
           rep =  substr(i, 102, 102),
           timestamp = as.POSIXct(TimeStamp, format = "%Y/%m/%d %H:%M:%S")) %>% 
    select(timestamp, shelf, rep, quantum = "Quantum..umol..m.2s..")
}))

data.time <- data %>% mutate(time = strftime(.$timestamp, format = "%H:%M:%S"),
                             date = strftime(.$timestamp, format = "%Y-%m-%d")) %>% 
  filter(time >= "10:00:00", time <= "13:00:00")


data.summary <- data.time %>% group_by(shelf) %>% 
  summarize(mean.quantum = mean(quantum))

ggplot(data.time, aes(x = factor(shelf), y = quantum)) + 
  geom_boxplot() +
  theme_bw()
