#### CLEAR THE ENVIRONMENT FIRST ---------------------------------------------

rm(list = ls(all.names = TRUE))


#### LOAD PACKAGES -----------------------------------------------------------

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


#### LOAD INCUBATION TEMPERATURE DATA ----------------------------------------

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
         block = factor(block))


#### FILTER TO EACH TRAITS' DATASET --------------------------------------------------------------

## filter to only eyed embryos
hatch.survival <- hatch %>% filter(eye != 0) %>% droplevels()

## filter to only hatched embryos
hatch.dpf <- hatch %>% filter(!is.na(dpf), hatch == 1) %>% droplevels()

## filter to only hatched embryos
hatch.ADD <- hatch %>% filter(!is.na(ADD), hatch == 1) %>% droplevels()


#### STATISTICAL ANALYSIS - SURVIVAL -----------------------------------------

## backward elimination to select best model
hatch.survival.glm <- buildmer(hatch ~ population + light + population:light
                                 (1|female:male) + (1|male) + (1|female) + (1|block) +
                                 (1|population:female) + (1|population:male) +
                                 (1|light:female) + (1|light:male), 
                               direction = 'backward', data = hatch.survival, 
                               family = binomial, control = glmerControl(optimizer = "bobyqa"))
( hatch.survival.glm.formula <- formula(hatch.survival.glm@model))

## fit best model
hatch.survival.glm.final <- glm(hatch.survival.glm.formula, data = hatch.survival, family = binomial)

## Calculate estimated marginal means - be very patient!
hatch.survival.glm.emm <- emmeans(hatch.survival.glm.final, ~ light | population)

## Pairwise
pairs(hatch.survival.glm.emm, simple = list("light", "population"), type = "response") 


#### STATISTICAL ANALYSIS - INCUBATION PERIOD (DPF) --------------------------

## fit full model
hatch.dpf.glm.full <- lmer(dpf ~ 1 + light + population + light:population + 
                             (1|female:male) + (1|male) + (1|female) + (1|block), 
                           data = hatch.dpf)

## backward elimination to select best model
hatch.dpf.glm <- step(hatch.dpf.glm.full)
( hatch.dpf.glm.formula <- get_model(hatch.dpf.glm)@call[["formula"]])

## fit best model
hatch.dpf.glm.final <- lmer(hatch.dpf.glm.formula, data = hatch.dpf)

## likelihood ratio test for fixed and random effects
mixed(hatch.dpf.glm.formula, data = hatch.dpf, method = "LRT")
rand(hatch.dpf.glm.final)


#### STATISTICAL ANALYSIS - INCUBATION PERIOD (ADD) --------------------------

## fit full model
hatch.ADD.glm.full <- lmer(ADD ~ 1 + light + population + light:population + 
                             (1|female:male) + (1|male) + (1|female) + (1|block), 
                           data = hatch.ADD)

## backward elimination to select best model
hatch.ADD.glm <- step(hatch.ADD.glm.full)
( hatch.ADD.glm.formula <- get_model(hatch.ADD.glm)@call[["formula"]])

## fit best model
hatch.ADD.glm.final <- lmer(hatch.ADD.glm.formula, data = hatch.ADD)

## likelihood ratio test for fixed and random effects
mixed(hatch.ADD.glm.formula, data = hatch.ADD, method = "LRT")
rand(hatch.ADD.glm.final)


#### CALCULATE MEAN AND SE FOR POPULATIONS -----------------------------------------------

## Embryo Survival Overall
hatch.survival.summary <- hatch %>% filter(eye != 0) %>% 
  group_by(population, light) %>% 
  summarize(mean.hatch = mean(hatch),
            se.hatch = sd(hatch)/sqrt(n()))

## Embryo Survival - Standardized Within Family
hatch.survival.summary.family <- hatch %>% filter(eye != 0) %>% 
  group_by(population, light, family) %>% 
  summarize(mean.hatch = mean(hatch)) %>% ungroup()

hatch.survival.stand <- hatch.survival.summary.family %>% filter(light == "Low") %>% 
  select(population, family, local.survival = mean.hatch)

hatch.survival.summary.stand <- hatch.survival.summary.family %>% left_join(hatch.survival.stand) %>% 
  mutate(survival.diff = 100*(1+(mean.hatch-local.survival)/local.survival)) %>%
  group_by(population, light) %>% 
  summarize(mean.survival.diff = mean(survival.diff),
            se.survival.diff = sd(survival.diff)/sqrt(n())) %>% 
  mutate(se.survival.diff = ifelse(se.survival.diff == 0, NA, se.survival.diff),
         percent.loss = 100-mean.survival.diff)


## Days Post Fertilization
hatch.dpf.summary <- hatch %>% filter(!is.na(dpf), hatch == 1) %>% 
  group_by(population, light) %>% 
  summarize(mean.dpf = mean(dpf),
            se.dpf = sd(dpf)/sqrt(n()))

## Days Post Fertilization - Standardized Within Family
hatch.dpf.summary.family <- hatch %>% filter(!is.na(dpf), hatch == 1) %>% 
  group_by(population, light, family) %>% 
  summarize(mean.dpf = mean(dpf)) %>% ungroup()

hatch.dpf.stand <- hatch.dpf.summary.family %>% filter(light == "Low") %>% 
  select(population, family, local.dpf = mean.dpf)

hatch.dpf.summary.stand <- hatch.dpf.summary.family %>% left_join(hatch.dpf.stand) %>% 
  mutate(dpf.diff = 100*(1+(mean.dpf-local.dpf)/local.dpf)) %>%
  group_by(population, light) %>% 
  summarize(mean.dpf.diff = mean(dpf.diff),
            se.dpf.diff = sd(dpf.diff)/sqrt(n())) %>% 
  mutate(se.dpf.diff = ifelse(se.dpf.diff == 0, NA, se.dpf.diff),
         percent.loss = 100-mean.dpf.diff)


## Accumulated Degree-Days
hatch.ADD.summary <- hatch %>% filter(!is.na(ADD), hatch == 1) %>% 
  group_by(population, light) %>% 
  summarize(mean.ADD = mean(ADD),
            se.ADD = sd(ADD)/sqrt(n()))

## Accumulated Degree-Days - Standardized Within Family
hatch.ADD.summary.family <- hatch %>% filter(!is.na(ADD), hatch == 1) %>% 
  group_by(population, light, family) %>% 
  summarize(mean.ADD = mean(ADD)) %>% ungroup()

hatch.ADD.stand <- hatch.ADD.summary.family %>% filter(light == "Low") %>% 
  select(population, family, local.ADD = mean.ADD)

hatch.ADD.summary.stand <- hatch.ADD.summary.family %>% left_join(hatch.ADD.stand) %>% 
  mutate(ADD.diff = 100*(1+(mean.ADD-local.ADD)/local.ADD)) %>%
  group_by(population, light) %>% 
  summarize(mean.ADD.diff = mean(ADD.diff),
            se.ADD.diff = sd(ADD.diff)/sqrt(n())) %>% 
  mutate(se.ADD.diff = ifelse(se.ADD.diff == 0, NA, se.ADD.diff),
         percent.loss = 100-mean.ADD.diff)

# VISUALIZATIONS - MEANS ----------------------------------------------------------------------

## Embryo Survival
plot.survival <- ggplot(hatch.survival.summary, aes(x = light, y = (mean.hatch * 100), group = population, color = population, shape = population, linetype = population)) + 
  geom_line(size = 1.0, position = position_dodge(0.13)) +
  geom_point(size = 5, position = position_dodge(0.13)) +
  geom_errorbar(aes(ymin = (mean.hatch - se.hatch) * 100, ymax = (mean.hatch + se.hatch) * 100), 
                position = position_dodge(0.13),
                size = 0.8, width = 0.1, linetype = "solid", show.legend = FALSE) +
  scale_x_discrete(labels = c("High", "Medium", "Low"), expand = c(0, 0.1)) +
  scale_y_continuous(limits = c(80, 101), breaks = seq(80, 100, 5), expand = c(0, 0)) +
  scale_color_grey("combine", start = 0.0, end = 0.6,
                   labels = c("Superior   ", "Ontario")) +
  scale_shape_manual("combine", values = c(1, 2), 
                     labels = c("Superior   ", "Ontario")) +
  scale_linetype_manual("combine", values = c("solid", "dashed"), 
                        labels = c("Superior   ", "Ontario")) +
  labs(y = "Mean Embryo Survival (%)") +
  theme_bw() +
  theme(axis.title.x = element_text(color = "Black", size = 22, margin = margin(10, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 22, margin = margin(0, 10, 0, 0)),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.ticks.length = unit(2, 'mm'),
        legend.title = element_blank(),
        legend.text = element_text(size = 20),
        legend.key.width = unit(1.25, 'cm'),
        legend.position = "top",
        plot.margin = unit(c(5, 5, 5, 5), 'mm'))

## Plot Standardized Survival
plot.survival.stand <- ggplot(hatch.survival.summary.stand, aes(x = population, y = mean.survival.diff, group = light, fill = light)) + 
  geom_bar(stat = "identity", size = 0.5, position = position_dodge(0.9), color = "black") +
  geom_errorbar(aes(ymin = (mean.survival.diff - se.survival.diff), ymax = (mean.survival.diff + se.survival.diff)), 
                position = position_dodge(0.9), size = 0.8, width = 0.4, show.legend = FALSE) +
  scale_y_continuous(limits = c(0.0, 130), breaks = seq(0.0, 130, 10), expand = c(0, 0)) +
  scale_fill_manual(values = c("#f4a582", "#92c5de", "#0571b0"),
                    labels = c("High  ", "Medium  ", "Low")) +
  coord_cartesian(ylim = c(80, 130)) +
  labs(y = "Standardized Survival (%)", x = "Population") +
  theme_bw() +
  theme(axis.title.x = element_text(color = "Black", size = 22, margin = margin(10, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 22, margin = margin(0, 10, 0, 0)),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.ticks.length = unit(2, 'mm'),
        legend.title = element_blank(),
        legend.text = element_text(size = 20),
        legend.key.size = unit(1.25, 'cm'),
        legend.position = "top",
        plot.margin = unit(c(5, 5, 5, 5), 'mm')) 

## Days Post Fertilization
plot.dpf <- ggplot(hatch.dpf.summary, aes(x = light, y = mean.dpf, group = population, color = population, shape = population, linetype = population)) + 
  geom_line(size = 1.0, position = position_dodge(0.13)) +
  geom_point(size = 5, position = position_dodge(0.13)) +
  geom_errorbar(aes(ymin = mean.dpf - se.dpf, ymax = mean.dpf + se.dpf), 
                position = position_dodge(0.13),
                size = 0.8, width = 0.1, linetype = "solid", show.legend = FALSE) +
  scale_x_discrete(labels = c("High", "Medium", "Low"), expand = c(0, 0.1)) +
  scale_y_continuous(limits = c(85, 130), breaks = seq(90, 130, 10), expand = c(0, 0)) +
  scale_color_grey("combine", start = 0.0, end = 0.6,
                   labels = c("Superior   ", "Ontario")) +
  scale_shape_manual("combine", values = c(1, 2), 
                     labels = c("Superior   ", "Ontario")) +
  scale_linetype_manual("combine", values = c("solid", "dashed"), 
                        labels = c("Superior   ", "Ontario")) +
  labs(y = "Mean DPF") +
  theme_bw() +
  theme(axis.title.x = element_text(color = "Black", size = 22, margin = margin(10, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 22, margin = margin(0, 10, 0, 0)),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.ticks.length = unit(2, 'mm'),
        legend.title = element_blank(),
        legend.text = element_text(size = 20),
        legend.key.width = unit(1.25, 'cm'),
        legend.position = "top",
        plot.margin = unit(c(5, 5, 5, 5), 'mm'))

## Plot Standardized DPF
plot.dpf.stand <- ggplot(hatch.dpf.summary.stand, aes(x = population, y = mean.dpf.diff, group = light, fill = light)) + 
  geom_bar(stat = "identity", size = 0.5, position = position_dodge(0.9), color = "black") +
  geom_errorbar(aes(ymin = (mean.dpf.diff - se.dpf.diff), ymax = (mean.dpf.diff + se.dpf.diff)), 
                position = position_dodge(0.9), size = 0.8, width = 0.4, show.legend = FALSE) +
  #scale_x_continuous(limits = c(1.75, 9.15), breaks = c(2, 4, 4.4, 6.9, 8, 8.9), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.0, 110), breaks = seq(0.0, 110, 10), expand = c(0, 0)) +
  scale_fill_manual(values = c("#f4a582", "#92c5de", "#0571b0"),
                    labels = c("High  ", "Medium  ", "Low")) +
  coord_cartesian(ylim = c(80, 110)) +
  labs(y = "Standardized DPF (%)", x = "Population") +
  theme_bw() +
  theme(axis.title.x = element_text(color = "Black", size = 22, margin = margin(10, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 22, margin = margin(0, 10, 0, 0)),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.ticks.length = unit(2, 'mm'),
        legend.title = element_blank(),
        legend.text = element_text(size = 20),
        legend.key.size = unit(1.25, 'cm'),
        legend.position = "top",
        plot.margin = unit(c(5, 5, 5, 5), 'mm')) 

## Accumulated Degree-Days
plot.ADD <- ggplot(hatch.ADD.summary, aes(x = light, y = mean.ADD, group = population, color = population, shape = population, linetype = population)) + 
  geom_line(size = 1.0, position = position_dodge(0.13)) +
  geom_point(size = 5, position = position_dodge(0.13)) +
  geom_errorbar(aes(ymin = mean.ADD - se.ADD, ymax = mean.ADD + se.ADD), 
                position = position_dodge(0.13),
                size = 0.8, width = 0.1, linetype = "solid", show.legend = FALSE) +
  scale_x_discrete(labels = c("High", "Medium", "Low"), expand = c(0, 0.1)) +
  scale_y_continuous(limits = c(350, 550), breaks = seq(350, 550, 50), expand = c(0, 0)) +
  scale_color_grey("combine", start = 0.0, end = 0.6,
                   labels = c("Superior   ", "Ontario")) +
  scale_shape_manual("combine", values = c(1, 2), 
                     labels = c("Superior   ", "Ontario")) +
  scale_linetype_manual("combine", values = c("solid", "dashed"), 
                        labels = c("Superior   ", "Ontario")) +
  labs(y = "Mean ADD (°C)") +
  theme_bw() +
  theme(axis.title.x = element_text(color = "Black", size = 22, margin = margin(10, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 22, margin = margin(0, 10, 0, 0)),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.ticks.length = unit(2, 'mm'),
        legend.title = element_blank(),
        legend.text = element_text(size = 20),
        legend.key.width = unit(1.25, 'cm'),
        legend.position = "top",
        plot.margin = unit(c(5, 5, 5, 5), 'mm'))

## Plot Standardized ADD
plot.ADD.stand <- ggplot(hatch.ADD.summary.stand, aes(x = population, y = mean.ADD.diff, group = light, fill = light)) + 
  geom_bar(stat = "identity", size = 0.5, position = position_dodge(0.9), color = "black") +
  geom_errorbar(aes(ymin = (mean.ADD.diff - se.ADD.diff), ymax = (mean.ADD.diff + se.ADD.diff)), 
                position = position_dodge(0.9), size = 0.8, width = 0.4, show.legend = FALSE) +
  scale_y_continuous(limits = c(0.0, 110), breaks = seq(0.0, 110, 10), expand = c(0, 0)) +
  scale_fill_manual(values = c("#f4a582", "#92c5de", "#0571b0"),
                    labels = c("High  ", "Medium  ", "Low")) +
  coord_cartesian(ylim = c(80, 110)) +
  labs(y = "Standardized ADD (%)", x = "Population") +
  theme_bw() +
  theme(axis.title.x = element_text(color = "Black", size = 22, margin = margin(10, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 22, margin = margin(0, 10, 0, 0)),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.ticks.length = unit(2, 'mm'),
        legend.title = element_blank(),
        legend.text = element_text(size = 20),
        legend.key.size = unit(1.25, 'cm'),
        legend.position = "top",
        plot.margin = unit(c(5, 5, 5, 5), 'mm')) 

## Combine all figures
plot.all <- grid.arrange(
  arrangeGrob(
    arrangeGrob(textGrob(""),
                get_legend(plot.survival),
                nrow = 1,
                widths = c(0.09, 1)),
    arrangeGrob(textGrob(""),
                get_legend(plot.survival.stand),
                nrow = 1,
                widths = c(0.09, 1)),
    ncol = 2,
    widths = c(1, 0.7)
  ),
  arrangeGrob(
    arrangeGrob(plot.survival + theme(legend.position = "none", axis.title.x = element_blank()),
                plot.dpf + theme(legend.position = "none", axis.title.x = element_blank()),
                plot.ADD + theme(legend.position = "none", axis.title.x = element_blank()),
                nrow = 3,
                bottom = textGrob("Light Treatment", x = 0.545, gp = gpar(cex = 2, fontfamily = "Arial"))),
    arrangeGrob(plot.survival.stand + theme(legend.position = "none", axis.title.x = element_blank()), 
                plot.dpf.stand + theme(legend.position = "none", axis.title.x = element_blank()),
                plot.ADD.stand + theme(legend.position = "none", axis.title.x = element_blank()),
                nrow = 3,
                bottom = textGrob("Population", x = 0.55, gp = gpar(cex = 2, fontfamily = "Arial"))),
    ncol = 2,
    widths = c(1, 0.7)
  ),
  heights = c(0.035, 1.1)
)

ggsave("figures/2020-Light-Embryo-LHT-SE.png", plot = plot.all, width = 18, height = 18, dpi = 200)

