#### LOAD PACKAGES -------------------------------------------------------------------------------

library(tidyverse)
library(readxl)
library(lme4)
library(lmerTest)
library(afex)
library(buildmer)
library(emmeans)
library(ggplot2)
library(gridExtra)
library(grid)
library(cowplot)


#### LOAD LARVAL LENGTH DATA ---------------------------------------------------------------------

larval.ls <- read_excel("data/Coregonine-Light-Experiment-LarvalMeasurements.xlsx", sheet = "Superior")
larval.lo <- read_excel("data/Coregonine-Light-Experiment-LarvalMeasurements.xlsx", sheet = "Ontario")

# Combine each population, temperature, and species
larval <- bind_rows(larval.ls, larval.lo) %>% 
  mutate(population = factor(population, levels = c("Superior", "Ontario"), ordered = TRUE),
         light = factor(light, ordered = TRUE, levels = c("Low", "Medium", "High")),
         female = factor(female, levels = seq(1, 12, 1),
                         labels = c("F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8", "F9", "F10", "F11", "F12")),
         male = factor(male, levels = seq(1, 16, 1),
                       labels = c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8", "M9", "M10", "M11", "M12", "M13", "M14", "M15", "M16")),
         block = factor(block),
         trans.yolk = y_vol_mm3^(1/3),
         trans.tl = length_mm^3)


#### FILTER TO EACH TRAITS' DATASET --------------------------------------------------------------

larval.tl <- larval %>% filter(!is.na(length_mm), length_mm != 0) %>% droplevels() %>% 
  filter(include.tl == "y")
larval.yolk <- larval %>% filter(!is.na(y_vol_mm3), y_vol_mm3 != 0) %>% droplevels() %>% 
  filter(include.yolk == "y")

## Clean up environment
#rm(larval.lo, larval.ls)


# STATISTICAL ANALYSIS - LENGTH-AT-HATCH - NA -------------------------------------------------

## fit full model
larval.tl.glm.full <- lmer(trans.tl ~ 1 + light + population + light:population + 
                             (1|family) + (1|male) + (1|female) + (1|block), 
                              data = larval.tl)

## backward elimination to select best model
larval.tl.glm <- step(larval.tl.glm.full)
( larval.tl.glm.formula <- get_model(larval.tl.glm)@call[["formula"]])

## fit best model
larval.tl.glm.final <- lmer(larval.tl.glm.formula, data = larval.tl)

## check residuals for normality
lattice::qqmath(larval.tl.glm.final, id = 0.1, idLabels = ~.obs)
hist(rstudent(larval.tl.glm.final))

## likelihood ratio test for fixed and random effects
mixed(larval.tl.glm.formula, data = larval.tl, method = "LRT")
rand(larval.tl.glm.final)

## Calculate estimated marginal means - be very patient!
larval.tl.glm.ontario <- lmer(trans.tl ~ light + (1 | male) + (1 | female), data = filter(larval.tl, population == "Ontario"))
larval.tl.glm.superior <- lmer(trans.tl ~ light + (1 | male) + (1 | female), data = filter(larval.tl, population == "Superior"))

larval.tl.glm.ontario.emm <- emmeans(larval.tl.glm.ontario, ~ light)
larval.tl.glm.superior.emm <- emmeans(larval.tl.glm.superior, ~ light)

## Pairwise
pairs(larval.tl.glm.ontario.emm, adjust = "tukey") 
pairs(larval.tl.glm.superior.emm, adjust = "tukey") 

#### STATISTICAL ANALYSIS - YOLK-SAC VOLUME ------------------------------------------------------

## fit full model
larval.yolk.glm.full <- lmer(trans.yolk ~ 1 + light + population + light:population + 
                                (1|family) + (1|male) + (1|female) + (1|block), 
                              data = larval.yolk)

## backward elimination to select best model
larval.yolk.glm <- step(larval.yolk.glm.full)
( larval.yolk.glm.formula <- get_model(larval.yolk.glm)@call[["formula"]])

## fit best model
larval.yolk.glm.final <- lmer(larval.yolk.glm.formula, data = larval.yolk)

## check residuals for normality
lattice::qqmath(larval.yolk.glm.final, id = 0.1, idLabels = ~.obs)
hist(rstudent(larval.yolk.glm.final))

## likelihood ratio test for fixed and random effects
mixed(larval.yolk.glm.formula, data = larval.yolk, method = "LRT")
rand(larval.yolk.glm.final)


#### CALCULATE MEAN AND SE FOR NA & FI POPULATIONS -----------------------------------------------

## Length-at-Hatch - Overall
larval.tl.summary <- larval.tl %>%
  group_by(population, light) %>% 
  summarize(mean.trait = mean(length_mm),
            se.trait = sd(length_mm)/sqrt(n())) %>% 
  mutate(trait = "LAH")

## Length-at-Hatch - Standardized Within Family
larval.tl.summary.family <- larval %>% filter(!is.na(length_mm), length_mm != 0) %>% 
  group_by(population, light, family) %>% 
  dplyr::summarize(mean.tl = mean(length_mm)) %>% ungroup()

larval.tl.stand <- larval.tl.summary.family %>% filter(light == "Low") %>% 
  select(population, family, local.tl = mean.tl)

larval.tl.summary.stand <- larval.tl.summary.family %>% left_join(larval.tl.stand) %>% 
  mutate(tl.diff = 100*((mean.tl-local.tl)/local.tl)) %>%
  filter(!is.na(tl.diff)) %>% 
  group_by(population, light) %>% 
  summarize(mean.trait.stand = mean(tl.diff),
            se.trait.stand = sd(tl.diff)/sqrt(n())) %>% 
  mutate(se.trait.stand = ifelse(se.trait.stand == 0, NA, se.trait.stand),
         trait = "LAH")
rm(larval.tl.summary.family, larval.tl.stand)

## Yolk-sac Volume - Overall
larval.yolk.summary <- larval.yolk %>% 
  group_by(population, light) %>% 
  summarize(mean.trait = mean(y_vol_mm3),
            se.trait = sd(y_vol_mm3)/sqrt(n())) %>% 
  mutate(trait = "YSV")

## Length-at-Hatch - Standardized Within Family
larval.yolk.summary.family <- larval %>% filter(!is.na(y_vol_mm3), y_vol_mm3 != 0) %>% 
  group_by(population, light, family) %>% 
  dplyr::summarize(mean.yolk = mean(y_vol_mm3)) %>% ungroup()

larval.yolk.stand <- larval.yolk.summary.family %>% filter(light == "Low") %>% 
  select(population, family, local.yolk = mean.yolk)

larval.yolk.summary.stand <- larval.yolk.summary.family %>% left_join(larval.yolk.stand) %>% 
  mutate(yolk.diff = 100*((mean.yolk-local.yolk)/local.yolk)) %>%
  filter(!is.na(yolk.diff)) %>% 
  group_by(population, light) %>% 
  summarize(mean.trait.stand = mean(yolk.diff),
            se.trait.stand = sd(yolk.diff)/sqrt(n())) %>% 
  mutate(se.trait.stand = ifelse(se.trait.stand == 0, NA, se.trait.stand),
         trait = "YSV")
rm(larval.yolk.summary.family, larval.yolk.stand)

