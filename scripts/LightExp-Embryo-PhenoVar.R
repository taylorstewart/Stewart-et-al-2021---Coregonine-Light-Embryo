#### LOAD PACKAGES & SET THREADS -----------------------------------------------------------------

library(tidyverse)
library(readxl)
library(data.table)
library(ggplot2)
library(fullfact)
library(parallel)
library(gridExtra)
library(grid)
library(cowplot)


#### LOAD INCUBATION TEMPERATURE DATA ------------------------------------------------------------

ADD <- read.csv("data/Artedi-Light-ADD-2020.csv", header = TRUE) %>% 
  dplyr::select(population, light, ADD) %>% 
  group_by(population, light) %>% 
  mutate(dpf = 1:n())


#### LOAD HATCHING DATA ------------------------------------------------------

hatch <- read_excel("data/Coregonine-Light-Experiment-Hatch.xlsx", sheet = "2020HatchingData") %>% 
  filter(include == "y", block != "A" | population != "superior") %>% 
  mutate(eye = as.numeric(eye),
         hatch = as.numeric(hatch)) %>% 
  left_join(ADD) %>% 
  dplyr::select(population, light, male, female, block, eye, hatch, dpf, ADD) %>% 
  mutate(population = factor(population, levels = c("Superior", "Ontario"), ordered = TRUE),
         light = factor(light, ordered = TRUE, levels = c("High", "Medium", "Low")),
         female = factor(female, levels = seq(1, 12, 1),
                         labels = c("F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8", "F9", "F10", "F11", "F12")),
         male = factor(male, levels = seq(1, 16, 1),
                       labels = c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8", "M9", "M10", "M11", "M12", "M13", "M14", "M15", "M16")),
         family = factor(paste0(female, male)),
         block = factor(block),
         group = interaction(population, light),
         group = gsub("°", "", group),
         group = gsub("\\.", "_", group)) %>% 
  rename(sire = male, dam = female)


#### FILTER TO EACH TRAITS' DATASET --------------------------------------------------------------

## filter to only eyed embryos
hatch.survival <- hatch %>% filter(eye != 0) %>% droplevels()

## filter to only hatched embryos
hatch.dpf <- hatch %>% filter(!is.na(dpf), hatch == 1) %>% droplevels()

## filter to only hatched embryos
hatch.ADD <- hatch %>% filter(!is.na(ADD), hatch == 1) %>% droplevels()


#### STATISTICAL ANALYSIS - GENERATE OBSERVED HERITABILITY ---------------------------------------

## Embryo Survival
phenoVar.survival.obs <- do.call(rbind, lapply(unique(hatch.survival$group), function(grp) {
  ## Filter to only a single group
  data.group <- hatch.survival %>% filter(group == grp) %>% 
    select(family, dam, sire, block, hatch)
  
  obs.survival <- observGlmer(observ = data.group, dam = "dam", sire = "sire", response = "hatch",
                              fam_link = binomial(logit))
  
  obs.survival.df <- data.frame(population = gsub("_", "", substr(grp, 1, 8)),
                                light = gsub("_", "", substr(grp, 9, nchar(grp))),
                                dam.var = obs.survival$random[3,2],
                                dam.p = obs.survival$random[3,7],
                                dam.perc = obs.survival$random[3,3],
                                sire.var = obs.survival$random[2,2],
                                sire.p = obs.survival$random[2,7],
                                sire.perc = obs.survival$random[2,3],
                                dam.sire.var = obs.survival$random[1,2],
                                dam.sire.p = obs.survival$random[1,7],
                                dam.sire.perc = obs.survival$random[1,3],
                                residual.var = obs.survival$other[1,2],
                                residual.perc = obs.survival$other[1,3]) %>% 
    mutate_if(is.numeric, round, 4)
})) %>% mutate(trait = "survival")

## DPF
phenoVar.dpf.obs <- do.call(rbind, lapply(unique(hatch.dpf$group), function(grp) {
  ## Filter to only a single group
  data.group <- hatch.dpf %>% filter(group == grp) %>% 
    select(family, dam, sire, block, dpf)
  
  obs.dpf <- observLmer(observ = data.group, dam = "dam", sire = "sire", response = "dpf")
  
  obs.dpf.df <- data.frame(population = gsub("_", "", substr(grp, 1, 8)),
                           light = gsub("_", "", substr(grp, 9, nchar(grp))),
                           dam.var = obs.dpf$random[3,2],
                           dam.p = obs.dpf$random[3,7],
                           dam.perc = obs.dpf$random[3,3],
                           sire.var = obs.dpf$random[2,2],
                           sire.p = obs.dpf$random[2,7],
                           sire.perc = obs.dpf$random[2,3],
                           dam.sire.var = obs.dpf$random[1,2],
                           dam.sire.p = obs.dpf$random[1,7],
                           dam.sire.perc = obs.dpf$random[1,3],
                           residual.var = obs.dpf$other[1,2],
                           residual.perc = obs.dpf$other[1,3]) %>% 
    mutate_if(is.numeric, round, 4)
})) %>% mutate(trait = "dpf")

## ADD
phenoVar.ADD.obs <- do.call(rbind, lapply(unique(hatch.ADD$group), function(grp) {
  ## Filter to only a single group
  data.group <- hatch.ADD %>% filter(group == grp) %>% 
    select(family, dam, sire, block, ADD)
  
  obs.ADD <- observLmer(observ = data.group, dam = "dam", sire = "sire", response = "ADD")
  
  obs.ADD.df <- data.frame(population = gsub("_", "", substr(grp, 1, 8)),
                           light = gsub("_", "", substr(grp, 9, nchar(grp))),
                           dam.var = obs.ADD$random[3,2],
                           dam.p = obs.ADD$random[3,7],
                           dam.perc = obs.ADD$random[3,3],
                           sire.var = obs.ADD$random[2,2],
                           sire.p = obs.ADD$random[2,7],
                           sire.perc = obs.ADD$random[2,3],
                           dam.sire.var = obs.ADD$random[1,2],
                           dam.sire.p = obs.ADD$random[1,7],
                           dam.sire.perc = obs.ADD$random[1,3],
                           residual.var = obs.ADD$other[1,2],
                           residual.perc = obs.ADD$other[1,3]) %>% 
    mutate_if(is.numeric, round, 4)  
})) %>% mutate(trait = "ADD")


#### CALCUALTE MEAN VARIANCE ACROSS TEMPERATURES -------------------------------------------------

phenoVar.embryo.mean <- bind_rows(phenoVar.survival.obs, phenoVar.dpf.obs, phenoVar.ADD.obs) %>% 
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

phenoVar.embryo.error <- bind_rows(phenoVar.survival.obs, phenoVar.dpf.obs, phenoVar.ADD.obs) %>% 
  group_by(population, trait) %>% 
  summarize(dam.perc.se = sd(dam.perc)/sqrt(n()),
            sire.perc.se = sd(sire.perc)/sqrt(n()),
            dam.sire.perc.se = sd(dam.sire.perc)/sqrt(n()),
            residual.perc.se = sd(residual.perc)/sqrt(n())) %>% 
  pivot_longer(3:6, names_to = "component", values_to = "error") %>% 
  mutate(component = factor(component, ordered = TRUE,
                            levels = c("dam.perc.se", "sire.perc.se", "dam.sire.perc.se", "residual.perc.se"),
                            labels = c("Dam", "Sire", "Dam.Sire", "Error")))


#### JOIN MEAN AND ERROR -------------------------------------------------------------------------

phenoVar.embryo.all <- left_join(phenoVar.embryo.mean, phenoVar.embryo.error) %>% 
  mutate(population = factor(population, ordered = TRUE, levels = c("Superior", "Ontario")))


#### CALCULATE CORRELATIONS ----------------------------------------------------------------------

## Embryo Survival
phenoVar.survival.cor <- phenoVar.survival.obs %>% 
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

## DPF
phenoVar.dpf.cor <- phenoVar.dpf.obs %>% 
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

## ADD
phenoVar.ADD.cor <- phenoVar.ADD.obs %>% 
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
