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
larval.tl.glm.emm <- emmeans(larval.tl.glm.final, ~ population)

## Pairwise
pairs(larval.tl.glm.emm, simple = list("population"), adjust = "tukey", type = "response") 


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
  summarize(mean.tl = mean(length_mm),
            se.tl = sd(length_mm)/sqrt(n()),
            n = n())

## Length-at-Hatch - Standardized Within Family
larval.tl.summary.family <- larval %>% filter(!is.na(length_mm), length_mm != 0) %>% 
  group_by(population, light, family) %>% 
  summarize(mean.tl = mean(length_mm)) %>% ungroup()

larval.tl.stand <- larval.tl.summary.family %>% filter(light == "Low") %>% 
  select(population, family, local.tl = mean.tl)

larval.tl.summary.stand <- larval.tl.summary.family %>% left_join(larval.tl.stand) %>% 
  mutate(tl.diff = 100*(1+(mean.tl-local.tl)/local.tl)) %>%
  filter(!is.na(tl.diff)) %>% 
  group_by(population, light) %>% 
  summarize(mean.tl.diff = mean(tl.diff),
            se.tl.diff = sd(tl.diff)/sqrt(n())) %>% 
  mutate(se.tl.diff = ifelse(se.tl.diff == 0, NA, se.tl.diff),
         percent.loss = 100-mean.tl.diff)


## Yolk-sac Volume - Overall
larval.yolk.summary <- larval.yolk %>% 
  group_by(population, light) %>% 
  summarize(mean.yolk = mean(y_vol_mm3),
            se.yolk = sd(y_vol_mm3)/sqrt(n()),
            n = n())

## Length-at-Hatch - Standardized Within Family
larval.yolk.summary.family <- larval %>% filter(!is.na(y_vol_mm3), y_vol_mm3 != 0) %>% 
  group_by(population, light, family) %>% 
  summarize(mean.yolk = mean(y_vol_mm3)) %>% ungroup()

larval.yolk.stand <- larval.yolk.summary.family %>% filter(light == "Low") %>% 
  select(population, family, local.yolk = mean.yolk)

larval.yolk.summary.stand <- larval.yolk.summary.family %>% left_join(larval.yolk.stand) %>% 
  mutate(yolk.diff = 100*(1+(mean.yolk-local.yolk)/local.yolk)) %>%
  filter(!is.na(yolk.diff)) %>% 
  group_by(population, light) %>% 
  summarize(mean.yolk.diff = mean(yolk.diff),
            se.yolk.diff = sd(yolk.diff)/sqrt(n())) %>% 
  mutate(se.yolk.diff = ifelse(se.yolk.diff == 0, NA, se.yolk.diff),
         percent.loss = 100-mean.yolk.diff)

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
  theme(axis.title.x = element_text(color = "Black", size = 22, margin = margin(10, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 22, margin = margin(0, 10, 0, 0)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.ticks.length = unit(2, 'mm'),
        legend.title = element_blank(),
        legend.text = element_text(size = 20),
        legend.key.size = unit(1.25, 'cm'),
        legend.position = "top",
        plot.margin = unit(c(5, 5, 5, 5), 'mm')) 

##
plot.tl.stand <- ggplot(larval.tl.summary.stand, aes(x = population, y = mean.tl.diff, group = light, fill = light)) + 
  geom_bar(stat = "identity", size = 0.5, position = position_dodge(0.9), color = "black") +
  geom_errorbar(aes(ymin = (mean.tl.diff - se.tl.diff), ymax = (mean.tl.diff + se.tl.diff)), 
                position = position_dodge(0.9), size = 0.8, width = 0.4, show.legend = FALSE) +
  #scale_x_continuous(limits = c(1.75, 9.15), breaks = c(2, 4, 4.4, 6.9, 8, 8.9), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.0, 105), breaks = seq(0.0, 105, 5), expand = c(0, 0)) +
  scale_fill_manual(values = c("#f4a582", "#92c5de", "#0571b0"),
                    labels = c("High  ", "Medium  ", "Low")) +
  coord_cartesian(ylim = c(80, 105)) +
  labs(y = "Standardized LAH (%)", x = "Population") +
  theme_bw() +
  theme(axis.title.x = element_text(color = "Black", size = 22, margin = margin(10, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 22, margin = margin(0, 10, 0, 0)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.ticks.length = unit(2, 'mm'),
        legend.title = element_blank(),
        legend.text = element_text(size = 20),
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
  theme(axis.title.x = element_text(color = "Black", size = 22, margin = margin(10, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 22, margin = margin(0, 10, 0, 0)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.ticks.length = unit(2, 'mm'),
        legend.title = element_blank(),
        legend.text = element_text(size = 20),
        legend.key.size = unit(1.25, 'cm'),
        legend.position = "top",
        plot.margin = unit(c(5, 5, 5, 5), 'mm')) 

## 
plot.yolk.stand <- ggplot(larval.yolk.summary.stand, aes(x = population, y = mean.yolk.diff, group = light, fill = light)) + 
  geom_bar(stat = "identity", size = 0.5, position = position_dodge(0.9), color = "black") +
  geom_errorbar(aes(ymin = (mean.yolk.diff - se.yolk.diff), ymax = (mean.yolk.diff + se.yolk.diff)), 
                position = position_dodge(0.9), size = 0.8, width = 0.4, show.legend = FALSE) +
  #scale_x_continuous(limits = c(1.75, 9.15), breaks = c(2, 4, 4.4, 6.9, 8, 8.9), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.0, 125), breaks = seq(0.0, 120, 10), expand = c(0, 0)) +
  scale_fill_manual(values = c("#f4a582", "#92c5de", "#0571b0"),
                    labels = c("High  ", "Medium  ", "Low")) +
  coord_cartesian(ylim = c(80, 125)) +
  labs(y = "Standardized YSV (%)", x = "Population") +
  theme_bw() +
  theme(axis.title.x = element_text(color = "Black", size = 22, margin = margin(10, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 22, margin = margin(0, 10, 0, 0)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
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
                get_legend(plot.tl),
                nrow = 1,
                widths = c(0.09, 1)),
    arrangeGrob(textGrob(""),
                get_legend(plot.tl.stand),
                nrow = 1,
                widths = c(0.09, 1)),
    ncol = 2,
    widths = c(1, 0.7)
  ),
  arrangeGrob(
    arrangeGrob(plot.tl + theme(legend.position = "none", axis.title.x = element_blank()),
                plot.yolk + theme(legend.position = "none", axis.title.x = element_blank()),
                nrow = 2,
                bottom = textGrob("Light Treatment", x = 0.545, gp = gpar(cex = 2, fontfamily = "Arial"))),
    arrangeGrob(plot.tl.stand + theme(legend.position = "none", axis.title.x = element_blank()), 
                plot.yolk.stand + theme(legend.position = "none", axis.title.x = element_blank()),
                nrow = 2,
                bottom = textGrob("Study Group", x = 0.55, gp = gpar(cex = 2, fontfamily = "Arial"))),
    ncol = 2,
    widths = c(1, 0.7)
  ),
  heights = c(0.04, 1.1)
)

ggsave("figures/2020-Light-Larval-MT-SE.png", plot = plot.all, width = 18, height = 14, dpi = 300)
