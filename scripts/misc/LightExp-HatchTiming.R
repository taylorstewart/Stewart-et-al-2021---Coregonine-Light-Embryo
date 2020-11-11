## ===========================================================
## Clear the environment first
## ===========================================================
rm(list = ls(all.names = TRUE))


## ===========================================================
## Load packages
## ===========================================================
library(tidyverse)
library(magrittr)
library(readxl)
library(lubridate)
library(lme4)


## ===========================================================
## Load and Initial Manipulations of the Larval Hatch 
##   Number of eggs hatched
## ===========================================================
hatch.data <- read_excel("data/Artedi-Light-Experiment-2020-v2.xlsx", sheet = "HatchingData") %>% 
  filter(is.na(notes) | notes != "empty well") %>% 
  filter(block != "1" | population != "superior", !is.na(dpf)) %>% 
  mutate(hatch = as.numeric(hatch),
         treatment = factor(treatment, levels = c("high", "medium", "low"), ordered = TRUE),
         dpf = as.numeric(dpf)) %>% 
  filter(premature == 0, hatch != 0) %>% 
  dplyr::select(population, family, male, female, treatment, hatch, dpf)


## ===========================================================
## Calculate cumulative total, total number, and proportion for each treatment
## ===========================================================
hatch.cumsum <- hatch.data %>% group_by(population, treatment, dpf) %>% 
  summarize(quantity = sum(hatch)) %>% 
  mutate(prop = cumsum(quantity)/sum(quantity)) %>% ungroup() %>% 
  mutate(group = factor(paste0(population, "-", treatment)))


## ===========================================================
## 
## ===========================================================
lme.superior <- glm(prop ~ 0 + dpf + treatment, data = filter(hatch.cumsum, population == "superior"), family = binomial)
lme.ontario <- glm(prop ~ 0 + dpf + treatment, data = filter(hatch.cumsum, population == "ontario"), family = binomial)
summary(lme.superior)
summary(lme.ontario)


(sup.H.50 <- as.Date(as.numeric((log(0.5 / (1 - 0.5)) - lme.superior$coefficients[2]) / lme.superior$coefficients[1]), origin = "2019-12-06"))
(sup.M.50 <- as.Date(as.numeric((log(0.5 / (1 - 0.5)) - lme.superior$coefficients[3]) / lme.superior$coefficients[1]), origin = "2019-12-06"))
(sup.L.50 <- as.Date(as.numeric((log(0.5 / (1 - 0.5)) - lme.superior$coefficients[4]) / lme.superior$coefficients[1]), origin = "2019-12-06"))

(ont.H.50 <- as.Date(as.numeric((log(0.5 / (1 - 0.5)) - lme.ontario$coefficients[2]) / lme.ontario$coefficients[1]), origin = "2019-11-26"))
(ont.M.50 <- as.Date(as.numeric((log(0.5 / (1 - 0.5)) - lme.ontario$coefficients[3]) / lme.ontario$coefficients[1]), origin = "2019-11-26"))
(ont.L.50 <- as.Date(as.numeric((log(0.5 / (1 - 0.5)) - lme.ontario$coefficients[4]) / lme.ontario$coefficients[1]), origin = "2019-11-26"))

hatch.50 <- data.frame(population = c(rep("ontario", 3), rep("superior", 3)),
                       temperature = rep(c("H", "M", "L"), 2),
                       hatch.50.date = c(ont.H.50, ont.M.50, ont.L.50,
                                         sup.H.50, sup.M.50, sup.L.50))

write.csv(hatch.50, "data/2020-Artedi-Hatch50.csv", row.names = FALSE)


## Logistic regression with the 50% hatched line highlighted
ggplot(hatch.cumsum, aes(x = dpf, y = prop, colour = group)) + 
  geom_point(show.legend = FALSE, size = 3.5) + 
  stat_smooth(aes(y = prop, group = group),  method = "glm", 
              method.args = list(family = "binomial"), se = FALSE, size = 2) +
  geom_hline(yintercept = 0.50, colour = "#4d4d4d", size = 0.65) + 
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  #scale_x_continuous(limits = c(225, 410), breaks = seq(225, 400, 25), expand = c(0, 0)) +
  #scale_color_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00",
  #                              "#fdbf6f", "#a65628", "#f781bf", "#999999")) +
  #labs(x = "Degree Days (°C)", y = "Proportion of Hatched Larvae", colour = "Female\nFamilies") +
  #annotate("text", x = 236.5, y = 0.95, label = "2.0°C", size = 10) +
  theme_bw() +
  theme(panel.background = element_blank(), panel.grid = element_blank(), 
        axis.line = element_line(), axis.ticks.length = unit(2, 'mm'),
        axis.text.y = element_text(size = 30, colour = "black"),
        axis.text.x = element_text(size = 30, colour = "black"),
        axis.title.y = element_text(size = 35, margin = margin(0, 20, 0, 0)),
        axis.title.x = element_text(size = 35, margin = margin(20, 0, 0, 0)),
        legend.title = element_text(size = 30), legend.text = element_text(size = 30), 
        legend.key.size = unit(1.3, 'cm'), legend.position = "top",
        plot.margin = unit(c(8, 5, 5, 5), "mm"))
