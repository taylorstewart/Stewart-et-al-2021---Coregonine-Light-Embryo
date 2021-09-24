#### CLEAR THE ENVIRONMENT FIRST -----------------------------------------------------------------

rm(list = ls(all.names = TRUE))


#### LOAD PACKAGES -------------------------------------------------------------------------------

library(ggplot2)
library(gridExtra)
library(grid)
library(cowplot)
library(egg)


# RUN EACH PHENO SCRIPT -----------------------------------------------------------------------

source("scripts/LightExp-Embryo-LHT.R")
source("scripts/LightExp-Larvae-MT.R")

traitsOverall.all <- bind_rows(hatch.survival.summary, hatch.dpf.summary, hatch.ADD.summary, larval.tl.summary, larval.yolk.summary)
traitsStand.all <- bind_rows(hatch.survival.summary.stand, hatch.dpf.summary.stand, hatch.ADD.summary.stand, larval.tl.summary.stand, larval.yolk.summary.stand)


#### VISUALIZATIONS - MEANS ----------------------------------------------------------------------

## Embryo Survival
plot.survival <- ggplot(filter(traitsOverall.all, trait == "survival"), aes(x = population, y = (mean.trait*100), group = light, color = light, shape = light)) + 
  geom_point(size = 5, position = position_dodge(0.6), stroke = 1.5) +
  geom_errorbar(aes(ymin = (mean.trait - se.trait) * 100, ymax = (mean.trait + se.trait) * 100), 
                position = position_dodge(0.6),
                size = 1, width = 0.25, linetype = "solid", show.legend = FALSE) +
  scale_y_continuous(limits = c(80, 101), breaks = seq(80, 100, 5), expand = c(0, 0)) +
  scale_shape_manual("combine", values = c(21, 22, 24), labels = c("Low ", "Medium ", "High")) +
  scale_color_manual("combine", values = c("#0571b0", "#92c5de", "#f4a582"),
                     labels = c("Low ", "Medium ", "High")) +
  labs(y = "Mean Embryo Survival (%)", x = "Population") +
  theme_bw() +
  theme(axis.title.x = element_text(color = "Black", size = 22, margin = margin(10, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 22, margin = margin(0, 10, 0, 0)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.ticks.length = unit(2, 'mm'),
        legend.title = element_blank(),
        legend.text = element_text(size = 20),
        legend.key.width = unit(1.25, 'cm'),
        legend.position = "top",
        plot.margin = unit(c(5, 5, 5, 5), 'mm'))

## Days Post Fertilization
plot.dpf <- ggplot(filter(traitsOverall.all, trait == "dpf"), aes(x = population, y = mean.trait, group = light, color = light, shape = light)) + 
  geom_point(size = 5, position = position_dodge(0.6), stroke = 1.5) +
  geom_errorbar(aes(ymin = mean.trait - se.trait, ymax = mean.trait + se.trait), 
                position = position_dodge(0.6),
                size = 1, width = 0.25, linetype = "solid", show.legend = FALSE) +
  scale_y_continuous(limits = c(95, 120), breaks = seq(95, 130, 5), expand = c(0, 0)) +
  scale_shape_manual("combine", values = c(21, 22, 24), labels = c("Low ", "Medium ", "High")) +
  scale_color_manual("combine", values = c("#0571b0", "#92c5de", "#f4a582"),
                     labels = c("Low ", "Medium ", "High")) +
  labs(y = "Mean Days Post-fertilization", x = "Population") +
  theme_bw() +
  theme(axis.title.x = element_text(color = "Black", size = 22, margin = margin(10, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 22, margin = margin(0, 10, 0, 0)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.ticks.length = unit(2, 'mm'),
        legend.title = element_blank(),
        legend.text = element_text(size = 20),
        legend.key.width = unit(1.25, 'cm'),
        legend.position = "top",
        plot.margin = unit(c(5, 5, 5, 5), 'mm'))

## Accumulated Degree-Days
plot.add <- ggplot(filter(traitsOverall.all, trait == "ADD"), aes(x = population, y = mean.trait, group = light, color = light, shape = light)) + 
  geom_point(size = 5, position = position_dodge(0.6), stroke = 1.5) +
  geom_errorbar(aes(ymin = mean.trait - se.trait, ymax = mean.trait + se.trait), 
                position = position_dodge(0.6),
                size = 1, width = 0.25, linetype = "solid", show.legend = FALSE) +
  scale_y_continuous(limits = c(400, 505), breaks = seq(350, 500, 25), expand = c(0, 0)) +
  scale_shape_manual("combine", values = c(21, 22, 24), labels = c("Low ", "Medium ", "High")) +
  scale_color_manual("combine", values = c("#0571b0", "#92c5de", "#f4a582"),
                     labels = c("Low ", "Medium ", "High")) +
  labs(y = "Mean Accumulated\nDegree-days (°C)", x = "Population") +
  theme_bw() +
  theme(axis.title.x = element_text(color = "Black", size = 22, margin = margin(10, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 22, margin = margin(0, 10, 0, 0)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.ticks.length = unit(2, 'mm'),
        legend.title = element_blank(),
        legend.text = element_text(size = 20),
        legend.key.width = unit(1.25, 'cm'),
        legend.position = "top",
        plot.margin = unit(c(5, 1, 5, 1), 'mm'))

## Length-at-Hatch
plot.tl <- ggplot(filter(traitsOverall.all, trait == "LAH"), aes(x = population, y = mean.trait, group = light, color = light, shape = light)) + 
  geom_point(size = 5, position = position_dodge(0.6), stroke = 1.5) +
  geom_errorbar(aes(ymin = mean.trait - se.trait, ymax = mean.trait + se.trait), 
                position = position_dodge(0.6),
                size = 1, width = 0.25, linetype = "solid", show.legend = FALSE) +
  scale_y_continuous(limits = c(9.5, 11), breaks = seq(9.5, 11, 0.5), expand = c(0, 0)) +
  scale_shape_manual("combine", values = c(21, 22, 24), labels = c("Low ", "Medium ", "High")) +
  scale_color_manual("combine", values = c("#0571b0", "#92c5de", "#f4a582"),
                     labels = c("Low ", "Medium ", "High")) +
  labs(y = "Mean Length-at-Hatch (mm)", x = "Population") +
  theme_bw() +
  theme(axis.title.x = element_text(color = "Black", size = 22, margin = margin(10, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 22, margin = margin(0, 10, 0, 0)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.ticks.length = unit(2, 'mm'),
        legend.title = element_blank(),
        legend.text = element_text(size = 20),
        legend.key.width = unit(1.25, 'cm'),
        legend.position = "top",
        plot.margin = unit(c(5, 5, 5, 5), 'mm'))

## Yolk-sac Volume
plot.ysv <- ggplot(filter(traitsOverall.all, trait == "YSV"), aes(x = population, y = mean.trait, group = light, color = light, shape = light)) + 
  geom_point(size = 5, position = position_dodge(0.6), stroke = 1.5) +
  geom_errorbar(aes(ymin = mean.trait - se.trait, ymax = mean.trait + se.trait), 
                position = position_dodge(0.6),
                size = 1, width = 0.25, linetype = "solid", show.legend = FALSE) +
  scale_y_continuous(limits = c(0.3, 0.7), breaks = seq(0.3, 0.7, 0.1), expand = c(0, 0)) +
  scale_shape_manual("combine", values = c(21, 22, 24), labels = c("Low ", "Medium ", "High")) +
  scale_color_manual("combine", values = c("#0571b0", "#92c5de", "#f4a582"),
                     labels = c("Low ", "Medium ", "High")) +
  labs(y = expression("Mean Yolk-sac Volume (mm"^3*")"), x = "Population") +
  theme_bw() +
  theme(axis.title.x = element_text(color = "Black", size = 22, margin = margin(10, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 22, margin = margin(0, 10, 0, 0)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.ticks.length = unit(2, 'mm'),
        legend.title = element_blank(),
        legend.text = element_text(size = 20),
        legend.key.width = unit(1.25, 'cm'),
        legend.position = "top",
        plot.margin = unit(c(5, 5, 5, 5), 'mm'))

## Combine all figures
plot.overall.all <- grid.arrange(
  arrangeGrob(get_legend(plot.survival)),
  arrangeGrob(plot.survival + theme(legend.position = "none", axis.title.x = element_blank()),
                plot.tl + theme(legend.position = "none", axis.title.x = element_blank()),
                plot.ysv + theme(legend.position = "none", axis.title.x = element_blank()),
                nrow = 1),
  arrangeGrob(textGrob("")),
  arrangeGrob(textGrob(""),
              plot.dpf + theme(legend.position = "none", axis.title.x = element_blank()), 
              plot.add + theme(legend.position = "none", axis.title.x = element_blank()),
              textGrob(""),
              nrow = 1,
              widths = c(0.5, 1, 1, 0.5)),
  arrangeGrob(textGrob("")),
  arrangeGrob(textGrob("Population", x = 0.5, just = "bottom", gp = gpar(cex = 2, fontfamily = "Arial"))),
  heights = c(0.15, 1.0, 0.05, 1.0, 0.05, 0.04)
)

ggsave("figures/summaryForDefense.tiff", plot = plot.overall.all, width = 14, height = 11, dpi = 250)


#### VISUALIZATIONS - STANDARDIZED ---------------------------------------------------------------

## Plot Standardized Survival
plot.survival.stand <- ggplot(filter(traitsStand.all, trait == "survival"), aes(x = population, y = mean.trait.stand, group = light, color = light, shape = light)) + 
  geom_point(size = 5, position = position_dodge(0.6), stroke = 1.5) +
  geom_hline(yintercept = 0, color = "gray30", linetype = "dashed", alpha = 0.5) +
  geom_errorbar(aes(ymin = (mean.trait.stand - se.trait.stand), ymax = (mean.trait.stand + se.trait.stand)), 
                position = position_dodge(0.6), size = 0.8, width = 0.4, show.legend = FALSE) +
  scale_y_continuous(limits = c(-20.0, 20), breaks = seq(-20.0, 20, 10), expand = c(0, 0)) +
  scale_shape_manual("combine", values = c(21, 22, 24), labels = c("Low ", "Medium ", "High")) +
  scale_color_manual("combine", values = c("#0571b0", "#92c5de", "#f4a582"),
                     labels = c("Low ", "Medium ", "High")) +
  labs(y = "Standardized Embryo Survival (%)", x = "Population") +
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

## Plot Standardized DPF
plot.dpf.stand <- ggplot(filter(traitsStand.all, trait == "dpf"), aes(x = population, y = mean.trait.stand, group = light, color = light, shape = light)) + 
  geom_point(size = 5, position = position_dodge(0.6), stroke = 1.5) +
  geom_hline(yintercept = 0, color = "gray30", linetype = "dashed", alpha = 0.5) +
  geom_errorbar(aes(ymin = (mean.trait.stand - se.trait.stand), ymax = (mean.trait.stand + se.trait.stand)), 
                position = position_dodge(0.6), size = 0.8, width = 0.4, show.legend = FALSE) +
  scale_y_continuous(limits = c(-20.0, 20), breaks = seq(-20.0, 20, 10), expand = c(0, 0)) +
  scale_shape_manual("combine", values = c(21, 22, 24), labels = c("Low ", "Medium ", "High")) +
  scale_color_manual("combine", values = c("#0571b0", "#92c5de", "#f4a582"),
                     labels = c("Low ", "Medium ", "High")) +
  labs(y = "Standardized Days\nPost-fertilization (%)", x = "Population") +
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
        plot.margin = unit(c(5, 1, 5, 1), 'mm')) 

## Plot Standardized ADD
plot.add.stand <- ggplot(filter(traitsStand.all, trait == "ADD"), aes(x = population, y = mean.trait.stand, group = light, color = light, shape = light)) + 
  geom_point(size = 5, position = position_dodge(0.6), stroke = 1.5) +
  geom_hline(yintercept = 0, color = "gray30", linetype = "dashed", alpha = 0.5) +
  geom_errorbar(aes(ymin = (mean.trait.stand - se.trait.stand), ymax = (mean.trait.stand + se.trait.stand)), 
                position = position_dodge(0.6), size = 0.8, width = 0.4, show.legend = FALSE) +
  scale_y_continuous(limits = c(-20.0, 20), breaks = seq(-20.0, 20, 10), expand = c(0, 0)) +
  scale_shape_manual("combine", values = c(21, 22, 24), labels = c("Low ", "Medium ", "High")) +
  scale_color_manual("combine", values = c("#0571b0", "#92c5de", "#f4a582"),
                     labels = c("Low ", "Medium ", "High")) +
  labs(y = "Standardized Accumulated\nDegree-days (%)", x = "Population") +
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
        plot.margin = unit(c(5, 1, 5, 1), 'mm')) 

## Plot Standardized LAH
plot.tl.stand <- ggplot(filter(traitsStand.all, trait == "LAH"), aes(x = population, y = mean.trait.stand, group = light, color = light, shape = light)) + 
  geom_point(size = 5, position = position_dodge(0.6), stroke = 1.5) +
  geom_hline(yintercept = 0, color = "gray30", linetype = "dashed", alpha = 0.5) +
  geom_errorbar(aes(ymin = (mean.trait.stand - se.trait.stand), ymax = (mean.trait.stand + se.trait.stand)), 
                position = position_dodge(0.6), size = 0.8, width = 0.4, show.legend = FALSE) +
  scale_y_continuous(limits = c(-20.0, 20), breaks = seq(-20.0, 20, 10), expand = c(0, 0)) +
  scale_shape_manual("combine", values = c(21, 22, 24), labels = c("Low ", "Medium ", "High")) +
  scale_color_manual("combine", values = c("#0571b0", "#92c5de", "#f4a582"),
                     labels = c("Low ", "Medium ", "High")) +
  labs(y = "Standardized\nLength-at-Hatch (%)", x = "Population") +
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
        plot.margin = unit(c(5, 1, 5, 1), 'mm')) 

## Plot Standardized YSV
plot.ysv.stand <- ggplot(filter(traitsStand.all, trait == "YSV"), aes(x = population, y = mean.trait.stand, group = light, color = light, shape = light)) + 
  geom_point(size = 5, position = position_dodge(0.6), stroke = 1.5) +
  geom_hline(yintercept = 0, color = "gray30", linetype = "dashed", alpha = 0.5) +
  geom_errorbar(aes(ymin = (mean.trait.stand - se.trait.stand), ymax = (mean.trait.stand + se.trait.stand)), 
                position = position_dodge(0.6), size = 0.8, width = 0.4, show.legend = FALSE) +
  scale_y_continuous(limits = c(-20.0, 20), breaks = seq(-20.0, 20, 10), expand = c(0, 0)) +
  scale_shape_manual("combine", values = c(21, 22, 24), labels = c("Low ", "Medium ", "High")) +
  scale_color_manual("combine", values = c("#0571b0", "#92c5de", "#f4a582"),
                     labels = c("Low ", "Medium ", "High")) +
  labs(y = "Standardized Yolk-sac\nVolumne (%)", x = "Population") +
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
        plot.margin = unit(c(5, 1, 5, 1), 'mm')) 

## Combine all figures
plot.stand.all <- grid.arrange(
  arrangeGrob(get_legend(plot.survival.stand)),
  arrangeGrob(plot.survival.stand + theme(legend.position = "none", axis.title.x = element_blank()),
              plot.tl.stand + theme(legend.position = "none", axis.title.x = element_blank()),
              plot.ysv.stand + theme(legend.position = "none", axis.title.x = element_blank()),
              nrow = 1),
  arrangeGrob(textGrob("")),
  arrangeGrob(textGrob(""),
              plot.dpf.stand + theme(legend.position = "none", axis.title.x = element_blank()), 
              plot.add.stand + theme(legend.position = "none", axis.title.x = element_blank()),
              textGrob(""),
              nrow = 1,
              widths = c(0.5, 1, 1, 0.5)),
  arrangeGrob(textGrob("")),
  arrangeGrob(textGrob("Population", x = 0.5, just = "bottom", gp = gpar(cex = 2, fontfamily = "Arial"))),
  heights = c(0.15, 1.0, 0.05, 1.0, 0.05, 0.04)
)

ggsave("figures/summaryForDefense-Stand.tiff", plot = plot.stand.all, width = 14, height = 11, dpi = 250)

