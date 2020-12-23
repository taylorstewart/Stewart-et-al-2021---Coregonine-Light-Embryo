#### LOAD PACKAGES -------------------------------------------------------------------------------

library(tidyverse)
library(readxl)
library(data.table)
library(ggplot2)
library(fullfact)
library(parallel)
library(gridExtra)
library(grid)
library(cowplot)


#### LOAD LARVAL LENGTH DATA ---------------------------------------------------------------------

larval.ls <- read_excel("data/Coregonine-Light-Experiment-LarvalMeasurements.xlsx", sheet = "Superior")
larval.lo <- read_excel("data/Coregonine-Light-Experiment-LarvalMeasurements.xlsx", sheet = "Ontario")

# Combine each population, temperature, and species
larval <- bind_rows(larval.ls, larval.lo) %>% 
  mutate(population = factor(population, levels = c("Superior", "Ontario"), ordered = TRUE),
         light = factor(light, ordered = TRUE, levels = c("High", "Medium", "Low")),
         female = factor(female, levels = seq(1, 12, 1),
                         labels = c("F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8", "F9", "F10", "F11", "F12")),
         male = factor(male, levels = seq(1, 16, 1),
                       labels = c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8", "M9", "M10", "M11", "M12", "M13", "M14", "M15", "M16")),
         block = factor(block),
         group = gsub("-", ".", interaction(population, light)),
         group = gsub("°", "", group),
         group = gsub("\\.", "_", group)) %>% 
  rename(sire = male, dam = female)

## Clean up environment
rm(larval.ls, larval.lo)


# FILTER TO EACH TRAITS' DATASET --------------------------------------------------------------

## filter out missing lengths
larval.tl <- larval %>% filter(!is.na(length_mm), length_mm != 0) %>% droplevels()

## filter out missing yolks
larval.yolk <- larval %>% filter(!is.na(y_vol_mm3), y_vol_mm3 != 0) %>% droplevels()


# STATISTICAL ANALYSIS - GENERATE OBSERVED HERITABILITY ---------------------------------------

## Length-at-Hatch
phenoVar.tl.obs <- do.call(rbind, lapply(unique(larval.tl$group), function(grp) {
  ## Filter to only a single temperature treatment
  data.grp <- larval.tl %>% filter(group == grp) %>% 
    select(family, dam, sire, block, length_mm)
  
  obs.tl <- observLmer(observ = data.grp, dam = "dam", sire = "sire", response = "length_mm")
  
  obs.tl.df <- data.frame(population = gsub("_", "", substr(grp, 1, 8)),
                          light = gsub("_", "", substr(grp, 9, nchar(grp))),
                          dam.var = obs.tl$random[3,2],
                          dam.p = obs.tl$random[3,7],
                          dam.perc = obs.tl$random[3,3],
                          sire.var = obs.tl$random[2,2],
                          sire.p = obs.tl$random[2,7],
                          sire.perc = obs.tl$random[2,3],
                          dam.sire.var = obs.tl$random[1,2],
                          dam.sire.p = obs.tl$random[1,7],
                          dam.sire.perc = obs.tl$random[1,3],
                          residual.var = obs.tl$other[1,2],
                          residual.perc = obs.tl$other[1,3]) %>% 
    mutate_if(is.numeric, round, 4)  
})) %>% mutate(trait = "LAH")

## Yolk-sac Volume
phenoVar.yolk.obs <- do.call(rbind, lapply(unique(larval.yolk$group), function(grp) {
  ## Filter to only a single temperature treatment
  data.grp <- larval.yolk %>% filter(group == grp) %>% 
    select(family, dam, sire, block, y_vol_mm3)
  
  obs.yolk <- observLmer2(observ = data.grp, dam = "dam", sire = "sire", response = "y_vol_mm3", block = "block")
  
  obs.yolk.df <- data.frame(population = gsub("_", "", substr(grp, 1, 8)),
                            light = gsub("_", "", substr(grp, 9, nchar(grp))),
                            dam.var = obs.yolk$random[3,2],
                            dam.p = obs.yolk$random[3,7],
                            dam.perc = obs.yolk$random[3,3],
                            sire.var = obs.yolk$random[2,2],
                            sire.p = obs.yolk$random[2,7],
                            sire.perc = obs.yolk$random[2,3],
                            dam.sire.var = obs.yolk$random[1,2],
                            dam.sire.p = obs.yolk$random[1,7],
                            dam.sire.perc = obs.yolk$random[1,3],
                            residual.var = obs.yolk$other[1,2],
                            residual.perc = obs.yolk$other[1,3]) %>% 
    mutate_if(is.numeric, round, 4)  
})) %>% mutate(trait = "YSV")


#### CALCUALTE MEAN VARIANCE ACROSS TEMPERATURES -------------------------------------------------

phenoVar.larval.mean <- bind_rows(phenoVar.tl.obs, phenoVar.yolk.obs) %>% 
  group_by(population, trait) %>% 
  summarize(dam.perc.mean = mean(dam.perc),
            sire.perc.mean = mean(sire.perc),
            dam.sire.perc.mean = mean(dam.sire.perc),
            residual.perc.mean = mean(residual.perc)) %>% 
  pivot_longer(3:6, names_to = "component", values_to = "variance") %>% 
  mutate(component = factor(component, ordered = TRUE,
                            levels = c("dam.perc.mean", "sire.perc.mean", "dam.sire.perc.mean", "residual.perc.mean"),
                            labels = c("Dam", "Sire", "Dam.Sire", "Error")))


#### CALCUALTE ERROR ACROSS TEMPERATURES ---------------------------------------------------------

phenoVar.larval.error <- bind_rows(phenoVar.tl.obs, phenoVar.yolk.obs) %>% 
  group_by(population, trait) %>% 
  summarize(dam.perc.se = sd(dam.perc)/sqrt(n()),
            sire.perc.se = sd(sire.perc)/sqrt(n()),
            dam.sire.perc.se = sd(dam.sire.perc)/sqrt(n()),
            residual.perc.se = sd(residual.perc)/sqrt(n())) %>% 
  pivot_longer(3:6, names_to = "component", values_to = "error") %>% 
  mutate(component = factor(component, ordered = TRUE,
                            levels = c("dam.perc.se", "sire.perc.se", "dam.sire.perc.se", "residual.perc.se"),
                            labels = c("Dam", "Sire", "Dam.Sire", "Error")))


# JOIN MEAN AND ERROR -------------------------------------------------------------------------

phenoVar.larval.all <- left_join(phenoVar.larval.mean, phenoVar.larval.error) %>% 
  mutate(population = factor(population, ordered = TRUE, levels = c("Superior", "Ontario")))


#### CALCULATE CORRELATIONS ----------------------------------------------------------------------

## LAH
phenoVar.tl.cor <- phenoVar.tl.obs %>% 
  group_by(population) %>% 
  summarize(dam.cor = cor(dam.perc, 1:3),
            sire.cor = cor(sire.perc, 1:3),
            dam.sire.cor = cor(dam.sire.perc, 1:3),
            error.cor = cor(residual.perc, 1:3)) %>% 
  mutate_if(is.numeric, round, 2) %>% 
  mutate(dam.cor.2 = ifelse(dam.cor >= 0.7, "POSITIVE", ifelse(dam.cor <= -0.7, "NEGATIVE", "NC")),
         sire.cor.2 = ifelse(sire.cor >= 0.7, "POSITIVE", ifelse(sire.cor <= -0.7, "NEGATIVE", "NC")),
         dam.sire.cor.2 = ifelse(dam.sire.cor >= 0.7, "POSITIVE", ifelse(dam.sire.cor <= -0.7, "NEGATIVE", "NC")),
         error.cor.2 = ifelse(error.cor >= 0.7, "POSITIVE", ifelse(error.cor <= -0.7, "NEGATIVE", "NC")))

## YSV
phenoVar.yolk.cor <- phenoVar.yolk.obs %>% 
  group_by(population) %>% 
  summarize(dam.cor = cor(dam.perc, 1:3),
            sire.cor = cor(sire.perc, 1:3),
            dam.sire.cor = cor(dam.sire.perc, 1:3),
            error.cor = cor(residual.perc, 1:3)) %>% 
  mutate_if(is.numeric, round, 2) %>% 
  mutate(dam.cor.2 = ifelse(dam.cor >= 0.7, "POSITIVE", ifelse(dam.cor <= -0.7, "NEGATIVE", "NC")),
         sire.cor.2 = ifelse(sire.cor >= 0.7, "POSITIVE", ifelse(sire.cor <= -0.7, "NEGATIVE", "NC")),
         dam.sire.cor.2 = ifelse(dam.sire.cor >= 0.7, "POSITIVE", ifelse(dam.sire.cor <= -0.7, "NEGATIVE", "NC")),
         error.cor.2 = ifelse(error.cor >= 0.7, "POSITIVE", ifelse(error.cor <= -0.7, "NEGATIVE", "NC")))
