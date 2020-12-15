#### CLEAR THE ENVIRONMENT FIRST -----------------------------------------------------------------

rm(list = ls(all.names = TRUE))


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


# COMBINE ALL TRAITS --------------------------------------------------------------------------

phenoVar.all <- bind_rows(phenoVar.tl.obs, phenoVar.yolk.obs) %>% 
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
         trait = factor(trait, ordered = TRUE, levels = c("LAH", "YSV"),
                        labels = c("Length-at-Hatch", "Yolk-sac Volume")),
         population = factor(population, ordered = TRUE, levels = c("Superior", "Ontario")))


#### CALCULATE CORRELATIONS ----------------------------------------------------------------------

## LAH
phenoVar.tl.cor <- phenoVar.tl.obs %>% 
  group_by(population) %>% 
  summarize(dam.cor = cor(dam.perc, temperature),
            sire.cor = cor(sire.perc, temperature),
            dam.sire.cor = cor(dam.sire.perc, temperature),
            error.cor = cor(residual.perc, temperature)) %>% 
  mutate_if(is.numeric, round, 2) %>% 
  mutate(dam.cor.2 = ifelse(dam.cor >= 0.7, "POSITIVE", ifelse(dam.cor <= -0.7, "NEGATIVE", "NC")),
         sire.cor.2 = ifelse(sire.cor >= 0.7, "POSITIVE", ifelse(sire.cor <= -0.7, "NEGATIVE", "NC")),
         dam.sire.cor.2 = ifelse(dam.sire.cor >= 0.7, "POSITIVE", ifelse(dam.sire.cor <= -0.7, "NEGATIVE", "NC")),
         error.cor.2 = ifelse(error.cor >= 0.7, "POSITIVE", ifelse(error.cor <= -0.7, "NEGATIVE", "NC")))

## YSV
phenoVar.yolk.cor <- phenoVar.yolk.obs %>% 
  group_by(population) %>% 
  summarize(dam.cor = cor(dam.perc, temperature),
            sire.cor = cor(sire.perc, temperature),
            dam.sire.cor = cor(dam.sire.perc, temperature),
            error.cor = cor(residual.perc, temperature)) %>% 
  mutate_if(is.numeric, round, 2) %>% 
  mutate(dam.cor.2 = ifelse(dam.cor >= 0.7, "POSITIVE", ifelse(dam.cor <= -0.7, "NEGATIVE", "NC")),
         sire.cor.2 = ifelse(sire.cor >= 0.7, "POSITIVE", ifelse(sire.cor <= -0.7, "NEGATIVE", "NC")),
         dam.sire.cor.2 = ifelse(dam.sire.cor >= 0.7, "POSITIVE", ifelse(dam.sire.cor <= -0.7, "NEGATIVE", "NC")),
         error.cor.2 = ifelse(error.cor >= 0.7, "POSITIVE", ifelse(error.cor <= -0.7, "NEGATIVE", "NC")))


#### VISUALIZATION - HERITABILITY --------------------------------------------

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

ggsave("figures/2020-Light-Larvae-PhenoVar.png", width = 10, height = 8.5, dpi = 300)
 
