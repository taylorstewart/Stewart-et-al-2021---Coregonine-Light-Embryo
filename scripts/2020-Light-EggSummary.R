## ===========================================================
## Clear the environment first
## ===========================================================
rm(list = ls(all.names = TRUE))


## ===========================================================
## Load packages
## ===========================================================
library(tidyverse)
library(readxl)
library(magrittr)
library(ggplot2)


## ===========================================================
## Load incubation temperature data
## ===========================================================
ADD <- read.csv("data/2020-Artedi-Light-ADD.csv", header = TRUE) %>% 
  dplyr::select(population, treatment, ADD) %>% 
  group_by(population, treatment) %>% 
  mutate(dpf = 1:n())


## ===========================================================
## Load hatching data
## ===========================================================
hatch <- read_excel("data/Artedi-Light-Experiment-2020-v2.xlsx", sheet = "HatchingData") %>% 
  filter(is.na(notes) | notes != "empty well") %>% 
  filter(block != "1" | population != "superior", !is.na(dpf)) %>% 
  #filter(treatment != "low" | population != "superior") %>%
  mutate(eye = as.numeric(eye),
         hatch = as.numeric(hatch),
         dpf = as.numeric(dpf)) %>% 
  filter(!is.na(eye), !is.na(hatch)) %>% 
  left_join(ADD) %>% 
  dplyr::select(population, family, male, female, treatment, eye, premature, hatch, dpf, ADD)


#----------------------------------------------------------------------------------------------#
################################    Results     #############################
#----------------------------------------------------------------------------------------------#

hatch.survival.summary <- hatch %>% filter(eye != 0) %>% group_by(population, treatment) %>% 
  summarize(mean.survival = (mean(hatch)*100),
            sd.survival = (sd(hatch)*100),
            se.survival = sd.survival / sqrt(n())) %>% 
  mutate(treatment = factor(treatment, levels = c("high", "medium", "low"), ordered = TRUE))

hatch.ADD.summary <- hatch %>% filter(!is.na(ADD), hatch == 1, premature != 1) %>% group_by(population, treatment) %>% 
  summarize(mean.ADD = mean(ADD),
            sd.ADD = sd(ADD),
            se.ADD = sd.ADD / sqrt(n())) %>% 
  mutate(treatment = factor(treatment, levels = c("high", "medium", "low"), ordered = TRUE))

hatch.dpf.summary <- hatch %>% filter(hatch == 1, premature != 1, !is.na(dpf)) %>% group_by(population, treatment) %>% 
  summarize(mean.dpf = mean(dpf),
            sd.dpf = sd(dpf),
            se.dpf = sd.dpf / sqrt(n())) %>% 
  mutate(treatment = factor(treatment, levels = c("high", "medium", "low"), ordered = TRUE))


ggplot(hatch.survival.summary, aes(x = treatment, y = mean.survival, group = population, color = population)) + 
  geom_line(size = 1.0, position = position_dodge(0.2)) +
  geom_point(size = 2.5, position = position_dodge(0.2)) +
  geom_errorbar(aes(ymin = mean.survival-se.survival, ymax = mean.survival+se.survival), size = 1.0, width = 0.2, position = position_dodge(0.2)) +
  scale_x_discrete(expand = c(0, 0.2), labels = c("High", "Medium", "Low")) +
  scale_y_continuous(limits = c(80, 100), breaks = seq(80, 100, 5), expand = c(0, 0)) +
  scale_color_manual(labels = c("L. Ontario    ", "L. Superior"),
                     values = c("#33a02c", "#1f78b4")) +
  labs(x = "Light Treatment (°C)", y = "Embryo Survival (% ± SE)", color = "Populations") +
  theme_classic() +
  theme(axis.title.x = element_text(color = "Black", size = 20, face = "bold", margin = margin(10, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 20, face = "bold", margin = margin(0, 10, 0, 0)),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        #legend.background = element_rect(size = 0.5, linetype = "solid", colour = "black"),
        #legend.title = element_text(face = "bold", size = 15),
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        legend.key.size = unit(1.0, 'cm'),
        legend.position = "top",
        #legend.position = c(0.17, 0.8),
        plot.margin = unit(c(5, 5, 5, 5), 'mm'))

ggsave("figures/2020-Light-Survival.png", width = 12, height = 7, dpi = 300)


ggplot(hatch.dpf.summary, aes(x = treatment, y = mean.dpf, group = population, color = population)) + 
  geom_line(size = 1.0, position = position_dodge(0.2)) +
  geom_point(size = 2.5, position = position_dodge(0.2)) +
  geom_errorbar(aes(ymin = mean.dpf-se.dpf, ymax = mean.dpf+se.dpf), size = 1.0, width = 0.2, position = position_dodge(0.2)) +
  scale_x_discrete(expand = c(0, 0.2), labels = c("High", "Medium", "Low")) +
  scale_y_continuous(limits = c(95, 120), breaks = seq(95, 120, 5), expand = c(0, 0)) +
  scale_color_manual(labels = c("L. Ontario    ", "L. Superior"),
                     values = c("#33a02c", "#1f78b4")) +
  labs(x = "Light Treatment (°C)", y = "Incubation Period (No. Days ± SE)", color = "Populations") +
  theme_classic() +
  theme(axis.title.x = element_text(color = "Black", size = 20, face = "bold", margin = margin(10, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 20, face = "bold", margin = margin(0, 10, 0, 0)),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        #legend.background = element_rect(size = 0.5, linetype = "solid", colour = "black"),
        #legend.title = element_text(face = "bold", size = 15),
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        legend.key.size = unit(1.0, 'cm'),
        legend.position = "top",
        #legend.position = c(0.8, 0.86),
        plot.margin = unit(c(5, 5, 5, 5), 'mm'))

ggsave("figures/2020-Light-DPF.png", width = 12, height = 7, dpi = 300)


ggplot(hatch.ADD.summary, aes(x = treatment, y = mean.ADD, group = population, color = population)) + 
  geom_line(size = 1.0, position = position_dodge(0.2)) +
  geom_point(size = 2.5, position = position_dodge(0.2)) +
  geom_errorbar(aes(ymin = mean.ADD-se.ADD, ymax = mean.ADD+se.ADD), size = 1.0, width = 0.2, position = position_dodge(0.2)) +
  scale_x_discrete(expand = c(0, 0.2), labels = c("High", "Medium", "Low")) +
  scale_y_continuous(limits = c(400, 500), breaks = seq(400, 500, 25), expand = c(0, 0)) +
  scale_color_manual(labels = c("L. Ontario    ", "L. Superior"),
                     values = c("#33a02c", "#1f78b4")) +
  labs(x = "Light Treatment (°C)", y = "Incubation Period (ADD °C ± SE)", color = "Populations") +
  theme_classic() +
  theme(axis.title.x = element_text(color = "Black", size = 20, face = "bold", margin = margin(10, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 20, face = "bold", margin = margin(0, 10, 0, 0)),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        #legend.background = element_rect(size = 0.5, linetype = "solid", colour = "black"),
        #legend.title = element_text(face = "bold", size = 15),
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        legend.key.size = unit(1.0, 'cm'),
        legend.position = "top",
        #legend.position = c(0.17, 0.88),
        plot.margin = unit(c(5, 5, 5, 5), 'mm'))

ggsave("figures/2020-Light-ADD.png", width = 12, height = 7, dpi = 300)


################################## statistique #########################

library(car)
library(lme4)
library(MuMIn)
library(ade4)
library(emmeans)

options(na.action = "na.fail")


################### survival #####################

hatch.survival <- hatch %>% filter(eye != 0)

c <- glmer(hatch ~ temperature + population + temperature * population +        # fixed
          (1|male) + (1|female),                                                # random
          family = binomial("logit"),
          data = hatch.survival,
          control = glmerControl(optimizer = "bobyqa"))  

dg1 <- dredge(c)                    # to select all model based on AICc
dg1

best <- get.models(dg1,"8")[[1]]    # select best model based on AICc
summary(best)

Anova(best)
Anova(best, type = "III")

# Post-hoc test:
best.emm <- emmeans(best, ~ temperature * population)
(best.emm.pair <- pairs(best.emm, simple = list("population", c("temperature")), adjust = "fdr"))
best.emm.pair.cld <- multcomp::cld(best.emm) %>% 
  select(population, temperature, group = .group)


##################### incubation  ########################

hatch.ADD <- hatch %>% filter(!is.na(ADD))

c <- lmer(ADD ~ population + temperature + population * temperature +
         (1|male) + (1|female), data = hatch.ADD, REML = FALSE)

dg1 <- dredge(c)                    # to select all model based on AICc
dg1

best <- get.models(dg1,"8")[[1]]    # select best model based on AICc
summary(best)

Anova(best)
Anova(best, type = "III")

# Post hoc test :
library(emmeans)
best.emm <- emmeans(best, ~ temperature * population)
pairs(best.emm, simple = list("population", c("temperature")), adjust = "fdr") 


