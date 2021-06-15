#### LOAD PACKAGES -----------------------------------------------------------

library(tidyverse)
library(readxl)
library(lme4)
library(lmerTest)
library(afex)
library(buildmer)
library(car)
library(emmeans)
library(ggplot2)
library(gridExtra)
library(grid)
library(cowplot)
library(egg)


#### LOAD INCUBATION TEMPERATURE DATA ----------------------------------------

ADD <- read.csv("data/Artedi-Light-ADD-2020.csv", header = TRUE) %>% 
  dplyr::select(population, light, ADD) %>% 
  group_by(population, light) %>% 
  mutate(dpf = 1:n())


#### LOAD HATCHING DATA ------------------------------------------------------

hatch <- read_excel("data/Coregonine-Light-Experiment-Hatch.xlsx", sheet = "2020HatchingData") %>% 
  filter(block != "A" | population != "superior") %>% 
  mutate(eye = as.numeric(eye),
         hatch = as.numeric(hatch)) %>% 
  left_join(ADD) %>% 
  dplyr::select(population, light, male, female, block, no, eye, hatch, dpf, ADD, include.survival, include.incubation) %>% 
  mutate(population = factor(population, levels = c("Superior", "Ontario"), ordered = TRUE),
         light = factor(light, ordered = TRUE, levels = c("Low", "Medium", "High")),
         female = factor(female, levels = seq(1, 12, 1),
                         labels = c("F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8", "F9", "F10", "F11", "F12")),
         male = factor(male, levels = seq(1, 16, 1),
                       labels = c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8", "M9", "M10", "M11", "M12", "M13", "M14", "M15", "M16")),
         family = factor(paste0(female, male)),
         block = factor(block),
         trans.dpf = dpf^(1/3),
         trans.ADD = ADD^(1/3))


#### FILTER TO EACH TRAITS' DATASET --------------------------------------------------------------

## filter to only eyed embryos
hatch.survival <- hatch %>% filter(eye != 0) %>% droplevels() %>% 
  filter(include.survival == "y")

## filter to only hatched embryos
hatch.dpf <- hatch %>% filter(!is.na(dpf), hatch == 1) %>% droplevels() %>% 
  filter(include.incubation == "y")

## filter to only hatched embryos
hatch.ADD <- hatch %>% filter(!is.na(ADD), hatch == 1) %>% droplevels() %>% 
  filter(include.incubation == "y")


#### STATISTICAL ANALYSIS - SURVIVAL -----------------------------------------

## backward elimination to select best model
hatch.survival.glm <- buildmer(hatch ~ population + light + population:light
                                 (1|female:male) + (1|male) + (1|female) + (1|block), 
                               direction = 'backward', data = hatch.survival, 
                               family = binomial, control = glmerControl(optimizer = "bobyqa"))
( hatch.survival.glm.formula <- formula(hatch.survival.glm@model))

## fit best model
hatch.survival.glm.final <- glm(hatch.survival.glm.formula, data = hatch.survival, family = binomial)

## likelihood ratio test for fixed effects
Anova(hatch.survival.glm.final, method = "LRT", type = "III")

## Calculate estimated marginal means - be very patient!
hatch.survival.glm.ontario <- glm(hatch ~ 1 + light, data = filter(hatch.survival, population == "Ontario"))
hatch.survival.glm.superior <- glm(hatch ~ 1 + light, data = filter(hatch.survival, population == "Superior"))

hatch.survival.glm.ontario.emm <- emmeans(hatch.survival.glm.ontario, ~ light)
hatch.survival.glm.superior.emm <- emmeans(hatch.survival.glm.superior, ~ light)

## Pairwise
pairs(hatch.survival.glm.ontario.emm, adjust = "tukey") 
pairs(hatch.survival.glm.superior.emm, adjust = "tukey") 


#### STATISTICAL ANALYSIS - INCUBATION PERIOD (DPF) --------------------------

## fit full model
hatch.dpf.glm.full <- lmer(trans.dpf ~ 1 + light + population + light:population + 
                             (1|female:male) + (1|male) + (1|female) + (1|block), 
                           data = hatch.dpf)

## backward elimination to select best model
hatch.dpf.glm <- step(hatch.dpf.glm.full)
( hatch.dpf.glm.formula <- get_model(hatch.dpf.glm)@call[["formula"]])

## fit best model
hatch.dpf.glm.final <- lmer(hatch.dpf.glm.formula, data = hatch.dpf)

## check residuals for normality
lattice::qqmath(hatch.dpf.glm.final, id = 0.1, idLabels = ~.obs)
hist(rstudent(hatch.dpf.glm.final))

## likelihood ratio test for fixed and random effects
mixed(hatch.dpf.glm.formula, data = hatch.dpf, method = "LRT")
rand(hatch.dpf.glm.final)


#### STATISTICAL ANALYSIS - INCUBATION PERIOD (ADD) --------------------------

## fit full model
hatch.ADD.glm.full <- lmer(trans.ADD ~ 1 + light + population + light:population + 
                             (1|female:male) + (1|male) + (1|female) + (1|block), 
                           data = hatch.ADD)

## backward elimination to select best model
hatch.ADD.glm <- step(hatch.ADD.glm.full)
( hatch.ADD.glm.formula <- get_model(hatch.ADD.glm)@call[["formula"]])

## fit best model
hatch.ADD.glm.final <- lmer(hatch.ADD.glm.formula, data = hatch.ADD)

## check residuals for normality
lattice::qqmath(hatch.ADD.glm.final, id = 0.1, idLabels = ~.obs)
hist(rstudent(hatch.ADD.glm.final))

## likelihood ratio test for fixed and random effects
mixed(hatch.ADD.glm.formula, data = hatch.ADD, method = "LRT")
rand(hatch.ADD.glm.final)


#### CALCULATE MEAN AND SE FOR POPULATIONS -----------------------------------------------

## Embryo Survival Overall
hatch.survival.summary <- hatch %>% filter(eye != 0) %>% 
  group_by(population, light) %>% 
  summarize(mean.trait = mean(hatch),
            se.trait = sd(hatch)/sqrt(n())) %>% 
  mutate(trait = "survival")

## Embryo Survival - Standardized Within Family
hatch.survival.summary.family <- hatch %>% filter(eye != 0) %>% 
  group_by(population, light, family) %>% 
  summarize(mean.hatch = mean(hatch)) %>% ungroup()

hatch.survival.stand <- hatch.survival.summary.family %>% filter(light == "Low") %>% 
  select(population, family, local.survival = mean.hatch)

hatch.survival.summary.stand <- hatch.survival.summary.family %>% left_join(hatch.survival.stand) %>% 
  mutate(survival.diff = 100*((mean.hatch-local.survival)/local.survival)) %>%
  group_by(population, light) %>% 
  summarize(mean.trait.stand = mean(survival.diff),
            se.trait.stand = sd(survival.diff)/sqrt(n())) %>% 
  mutate(se.trait.stand = ifelse(se.trait.stand == 0, NA, se.trait.stand),
         trait = "survival")
rm(hatch.survival.summary.family, hatch.survival.stand)


## Days Post Fertilization
hatch.dpf.summary <- hatch %>% filter(!is.na(dpf), hatch == 1) %>% 
  group_by(population, light) %>% 
  summarize(mean.trait = mean(dpf),
            se.trait = sd(dpf)/sqrt(n())) %>% 
  mutate(trait = "dpf")

## Days Post Fertilization - Standardized Within Family
hatch.dpf.summary.family <- hatch %>% filter(!is.na(dpf), hatch == 1) %>% 
  group_by(population, light, family) %>% 
  summarize(mean.dpf = mean(dpf)) %>% ungroup()

hatch.dpf.stand <- hatch.dpf.summary.family %>% filter(light == "Low") %>% 
  select(population, family, local.dpf = mean.dpf)

hatch.dpf.summary.stand <- hatch.dpf.summary.family %>% left_join(hatch.dpf.stand) %>% 
  mutate(dpf.diff = 100*((mean.dpf-local.dpf)/local.dpf)) %>%
  group_by(population, light) %>% 
  summarize(mean.trait.stand = mean(dpf.diff),
            se.trait.stand = sd(dpf.diff)/sqrt(n())) %>% 
  mutate(se.trait.stand = ifelse(se.trait.stand == 0, NA, se.trait.stand),
         trait = "dpf")
rm(hatch.dpf.summary.family, hatch.dpf.stand)

## Accumulated Degree-Days
hatch.ADD.summary <- hatch %>% filter(!is.na(ADD), hatch == 1) %>% 
  group_by(population, light) %>% 
  summarize(mean.trait = mean(ADD),
            se.trait = sd(ADD)/sqrt(n())) %>% 
  mutate(trait = "ADD")

## Accumulated Degree-Days - Standardized Within Family
hatch.ADD.summary.family <- hatch %>% filter(!is.na(ADD), hatch == 1) %>% 
  group_by(population, light, family) %>% 
  summarize(mean.ADD = mean(ADD)) %>% ungroup()

hatch.ADD.stand <- hatch.ADD.summary.family %>% filter(light == "Low") %>% 
  select(population, family, local.ADD = mean.ADD)

hatch.ADD.summary.stand <- hatch.ADD.summary.family %>% left_join(hatch.ADD.stand) %>% 
  mutate(ADD.diff = 100*((mean.ADD-local.ADD)/local.ADD)) %>%
  group_by(population, light) %>% 
  summarize(mean.trait.stand = mean(ADD.diff),
            se.trait.stand = sd(ADD.diff)/sqrt(n())) %>% 
  mutate(se.trait.stand = ifelse(se.trait.stand == 0, NA, se.trait.stand),
         trait = "ADD")
rm(hatch.ADD.summary.family, hatch.ADD.stand)

