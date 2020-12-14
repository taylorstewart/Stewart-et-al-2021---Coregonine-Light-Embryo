#### CLEAR THE ENVIRONMENT FIRST -----------------------------------------------------------------

rm(list = ls(all.names = TRUE))


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


# COMBINE ALL TRAITS --------------------------------------------------------------------------

phenoVar.all <- bind_rows(phenoVar.survival.obs, phenoVar.dpf.obs, phenoVar.ADD.obs) %>% 
  select(population, light, trait, dam.perc, sire.perc, dam.sire.perc, residual.perc) %>% 
  pivot_longer(4:7, names_to = "component", values_to = "variance") %>% 
  mutate(component = factor(component, ordered = TRUE,
                            levels = c("dam.perc", "sire.perc", "dam.sire.perc", "residual.perc"),
                            labels = c("Dam", "Sire", "Dam.Sire", "Error")),
         component.trt = factor(interaction(component, light), ordered = TRUE,
                                levels = c("Dam.High", "Dam.Medium", "Dam.Low",
                                           "Sire.High", "Sire.Medium", "Sire.Low",
                                           "Dam.Sire.High", "Dam.Sire.Medium", "Dam.Sire.Low",
                                           "Error.High", "Error.Medium", "Error.Low")),
         trait = factor(trait, ordered = TRUE, levels = c("survival", "dpf", "ADD"),
                        labels = c("Embryo Survival", "Incubation Period (DPF)", "Incubation Period (ADD)")),
         population = factor(population, ordered = TRUE, levels = c("Superior", "Ontario")))


#### VISUALIZATION - HERITABILITY ----------------------------------------------------------------

ggplot(phenoVar.all, aes(x = population, y = variance, group = component.trt, fill = component)) + 
  geom_bar(stat = "identity", size = 0.5, position = position_dodge(0.9), color = "black") +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0.5)) +
  scale_fill_manual(values = c("#7bccc4", "#f0f9e8", "#bae4bc", "#2b8cbe"),
                    labels = c("Dam  ", "Sire  ", "Dam x Sire  ", "Error")) +
  annotation_custom(textGrob("High-Medium-Low", gp = gpar(fontsize = 15, col = "grey30")), 
                    xmin = 1.5, xmax = 1.5, ymin = -7.5, ymax = -7.5) +
  coord_cartesian(clip = "off") +
  labs(y = "% of Total Phenotypic Variation", x = "Population\nLight Treatment") +
  theme_bw() +
  theme(axis.title.x = element_text(color = "Black", size = 22, margin = margin(10, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 22, margin = margin(0, 10, 0, 0)),
        axis.text.x = element_text(size = 16, margin = margin(2, 0, 25, 0)),
        axis.text.y = element_text(size = 16),
        axis.ticks.length = unit(2, 'mm'),
        legend.title = element_blank(),
        legend.text = element_text(size = 20),
        legend.key.size = unit(1.25, 'cm'),
        legend.position = "top",
        strip.text = element_text(size = 16),
        strip.background = element_rect(color = "white", fill = "white"),
        plot.margin = unit(c(1, 1, 1, 1), 'mm')) +
  facet_wrap(~trait)

ggsave("figures/2020-Light-Embryo-PhenoVar.png", width = 14, height = 8.5, dpi = 300)
