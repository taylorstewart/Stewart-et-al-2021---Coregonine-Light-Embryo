## ===========================================================
## Clear the environment first
## ===========================================================
rm(list = ls(all.names = TRUE))


## ===========================================================
## Load packages
## ===========================================================
library(tidyverse)
library(readxl)
library(lubridate)
library(magrittr)
library(ggplot2)
library(car)
library(lme4)
library(MuMIn)
library(ade4)
library(emmeans)

options(na.action = "na.fail")


## ===========================================================
## Load and Initial Manipulations of the Larval Hatch Length Data
## ===========================================================
larvae.tl.lo.h <- read_excel("data/Artedi-Light-Experiments-Length-2020.xlsx", sheet = "LO-H")[1:20]
larvae.tl.lo.m <- read_excel("data/Artedi-Light-Experiments-Length-2020.xlsx", sheet = "LO-M")[1:20]
larvae.tl.lo.l <- read_excel("data/Artedi-Light-Experiments-Length-2020.xlsx", sheet = "LO-L")[1:20]
larvae.tl.ls.h <- read_excel("data/Artedi-Light-Experiments-Length-2020.xlsx", sheet = "LS-H")[1:20]
larvae.tl.ls.m <- read_excel("data/Artedi-Light-Experiments-Length-2020.xlsx", sheet = "LS-M")[1:20]
larvae.tl.ls.l <- read_excel("data/Artedi-Light-Experiments-Length-2020.xlsx", sheet = "LS-L")[1:20]

larvae.tl <- bind_rows(larvae.tl.lo.h, larvae.tl.lo.m, larvae.tl.lo.l, 
                       larvae.tl.ls.h, larvae.tl.ls.m, larvae.tl.ls.l) %>% 
  filter(!is.na(length_mm)) %>% 
  mutate(population = factor(population),
         treatment = factor(treatment))

## -----------------------------------------------------------
## Calculate means and std. error for each treatment
## -----------------------------------------------------------
larvae.tl.summary <- larvae.tl %>% group_by(population, treatment) %>% 
  summarize(mean.tl = mean(length_mm),
            sd.tl = sd(length_mm),
            n = n(),
            se.tl = sd.tl/sqrt(n)) %>% 
  arrange(population, treatment) %>% 
  mutate(treatment = factor(treatment, ordered = TRUE, levels = c("high", "medium", "low")))


ggplot(larvae.tl.summary, aes(x = treatment, y = mean.tl, group = population, color = population)) + 
  geom_point(size = 3) + 
  geom_line(size = 1) +
  geom_errorbar(aes(ymin = mean.tl - se.tl, ymax = mean.tl + se.tl), width = 0.09, size = 0.85) +
  scale_y_continuous(limits = c(9.0, 11), breaks = seq(9.0, 11.0, 0.5), expand = c(0, 0)) +
  scale_x_discrete(expand = c(0.0, 0.1), labels = c("High", "Medium", "Low")) +
  scale_color_manual(labels = c("L. Ontario    ", "L. Superior"), 
                    values = c("#fc8d59", "#91bfdb")) +  
  labs(y = "Mean Length-at-Hatch (mm ± SE)", x = 'Photoperiod Intensity') +
  theme_classic() +
  theme(axis.title.x = element_text(color = "Black", size = 20, face = "bold", margin = margin(10, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 20, face = "bold", margin = margin(0, 10, 0, 0)),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.ticks.length = unit(1.5, "mm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        legend.key.size = unit(1.0, 'cm'),
        legend.position = "top",
        plot.margin = unit(c(5, 5, 5, 5), 'mm'))

ggsave("figures/2020-Light-LAH.png", width = 12, height = 7, dpi = 300)


##############################################################
## ANALYSIS
##############################################################
## -----------------------------------------------------------
## Fit model
## -----------------------------------------------------------
glm <- lmer(length_mm ~ population + treatment + population * treatment + (1|male) + (1|female), 
            data = larvae.tl, REML = FALSE)

dg1 <- dredge(glm)                    # to select all model based on AICc
dg1

best <- get.models(dg1,"8")[[1]]    # select best model based on AICc
summary(best)

Anova(best)
Anova(best, type = "III")

# Post-hoc test:
best.emm <- emmeans(best, ~ treatment * population)
(best.emm.pair <- pairs(best.emm, simple = list("population", c("treatment")), adjust = "fdr"))
best.emm.pair.cld <- multcomp::cld(best.emm) %>% 
  select(population, treatment, group = .group)

