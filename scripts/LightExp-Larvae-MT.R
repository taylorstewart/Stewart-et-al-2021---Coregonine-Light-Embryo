#### CLEAR THE ENVIRONMENT FIRST -----------------------------------------------------------------

rm(list = ls(all.names = TRUE))


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
         light = factor(light, ordered = TRUE, levels = c("High", "Medium", "Low")),
         female = factor(female, levels = seq(1, 12, 1),
                         labels = c("F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8", "F9", "F10", "F11", "F12")),
         male = factor(male, levels = seq(1, 16, 1),
                       labels = c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8", "M9", "M10", "M11", "M12", "M13", "M14", "M15", "M16")),
         block = factor(block))


#### FILTER TO EACH TRAITS' DATASET --------------------------------------------------------------

larval.tl <- larval %>% filter(!is.na(length_mm), length_mm != 0) %>% droplevels()
larval.yolk <- larval %>% filter(!is.na(y_vol_mm3), y_vol_mm3 != 0) %>% droplevels()

## Clean up environment
rm(larval.lo, larval.ls)


# STATISTICAL ANALYSIS - LENGTH-AT-HATCH - NA -------------------------------------------------

## fit full model
larval.tl.glm.full <- lmer(length_mm ~ 1 + light + population + light:population + 
                             (1|family) + (1|male) + (1|female) + (1|block), 
                              data = larval.tl)

## backward elimination to select best model
larval.tl.glm <- step(larval.tl.glm.full)
( larval.tl.glm.formula <- get_model(larval.tl.glm)@call[["formula"]])

## fit best model
larval.tl.glm.final <- lmer(larval.tl.glm.formula, data = larval.tl)

## likelihood ratio test for fixed and random effects
mixed(larval.tl.glm.formula, data = larval.tl, method = "LRT")
rand(larval.tl.glm.final)

## Calculate estimated marginal means - be very patient!
larval.tl.glm.emm <- emmeans(larval.tl.glm.final, ~ population)

## Pairwise
pairs(larval.tl.glm.emm, simple = list("population"), adjust = "tukey", type = "response") 


#### STATISTICAL ANALYSIS - YOLK-SAC VOLUME ------------------------------------------------------

## fit full model
larval.yolk.glm.full <- lmer(y_vol_mm3 ~ 1 + light + population + light:population + 
                                (1|family) + (1|male) + (1|female) + (1|block), 
                              data = larval.yolk)

## backward elimination to select best model
larval.yolk.glm <- step(larval.yolk.glm.full)
( larval.yolk.glm.formula <- get_model(larval.yolk.glm)@call[["formula"]])

## fit best model
larval.yolk.glm.final <- lmer(larval.yolk.glm.formula, data = larval.yolk)

## likelihood ratio test for fixed and random effects
mixed(larval.yolk.glm.formula, data = larval.yolk, method = "LRT")
rand(larval.yolk.glm.final)


#### CALCULATE MEAN AND SE FOR NA & FI POPULATIONS -----------------------------------------------

## Length-at-Hatch
larval.tl.summary <- larval.tl %>% 
  group_by(population, light) %>% 
  summarize(mean.tl = mean(length_mm),
            se.tl = sd(length_mm)/sqrt(n()),
            n = n())

## Yolk-sac Volume
larval.yolk.summary <- larval.yolk %>% 
  group_by(population, light) %>% 
  summarize(mean.yolk = mean(y_vol_mm3),
            se.yolk = sd(y_vol_mm3)/sqrt(n()),
            n = n())


#### VISUALIZATIONS ----------------------------------------------------------

## Length-at-Hatch
plot.tl <- ggplot(larval.tl.summary, aes(x = light, y = mean.tl, group = population, color = population, shape = population, linetype = population)) + 
  geom_line(size = 1.0, position = position_dodge(0)) +
  geom_point(size = 3.25, position = position_dodge(0)) +
  geom_errorbar(aes(ymin = mean.tl - se.tl, ymax = mean.tl + se.tl), 
                position = position_dodge(0),
                size = 0.8, width = 0.1, linetype = "solid", show.legend = FALSE) +
  scale_x_discrete(labels = c("High", "Medium", "Low"), expand = c(0, 0.1)) +
  scale_y_continuous(limits = c(9.4, 11), breaks = seq(9.5, 11, 0.5), expand = c(0, 0)) +
  scale_color_grey("combine", start = 0.0, end = 0.6,
                   labels = c("Superior   ", "Ontario")) +
  scale_shape_manual("combine", values = c(1, 2), 
                     labels = c("Superior   ", "Ontario")) +
  scale_linetype_manual("combine", values = c("solid", "dashed"), 
                        labels = c("Superior   ", "Ontario")) +
  labs(y = "Mean LAH (mm)", x = "Light Treatment") +
  theme_bw() +
  theme(axis.title.x = element_text(color = "Black", size = 18, margin = margin(10, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 18, margin = margin(0, 10, 0, 0)),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        legend.key.size = unit(1.25, 'cm'),
        legend.position = "top",
        plot.margin = unit(c(5, 5, 5, 5), 'mm'))

## Yolk-sac Volume
plot.yolk <- ggplot(larval.yolk.summary, aes(x = light, y = mean.yolk, group = population, color = population, shape = population, linetype = population)) + 
  geom_line(size = 1.0, position = position_dodge(0)) +
  geom_point(size = 3.25, position = position_dodge(0)) +
  geom_errorbar(aes(ymin = mean.yolk - se.yolk, ymax = mean.yolk + se.yolk), 
                position = position_dodge(0),
                size = 0.8, width = 0.1, linetype = "solid", show.legend = FALSE) +
  scale_x_discrete(labels = c("High", "Medium", "Low"), expand = c(0, 0.1)) +
  scale_y_continuous(limits = c(0.3, 0.7), breaks = seq(0.3, 0.7, 0.1), expand = c(0, 0)) +
  scale_color_grey("combine", start = 0.0, end = 0.6,
                   labels = c("Superior   ", "Ontario")) +
  scale_shape_manual("combine", values = c(1, 2), 
                     labels = c("Superior   ", "Ontario")) +
  scale_linetype_manual("combine", values = c("solid", "dashed"), 
                        labels = c("Superior   ", "Ontario")) +
  labs(y = expression("Mean YSV (mm"^3*")"), x = "Light Treatment") +
  theme_bw() +
  theme(axis.title.x = element_text(color = "Black", size = 18, margin = margin(10, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 18, margin = margin(0, 10, 0, 0)),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        legend.key.size = unit(1.25, 'cm'),
        legend.position = "top",
        plot.margin = unit(c(5, 5, 5, 5), 'mm'))

## Combine all figures
plot.all <- grid.arrange(arrangeGrob(textGrob(""), 
                                     get_legend(plot.tl),
                                     nrow = 1,
                                     widths = c(0.03, 1)),
                         arrangeGrob(plot.tl + theme(legend.position = "none", axis.title.x = element_blank()), 
                                     plot.yolk + theme(legend.position = "none", axis.title.x = element_blank()),
                                     nrow = 2,
                                     bottom = textGrob("Light Treatment", x = 0.545, gp = gpar(cex = 1.5, fontfamily = "Arial"))),
                         heights = c(0.025, 1)
                         )

ggsave("figures/2020-Light-Larval-MT-SE.png", plot = plot.all, width = 8, height = 10, dpi = 300)

