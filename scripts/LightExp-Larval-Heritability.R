#### CLEAR THE ENVIRONMENT FIRST -----------------------------------------------------------------

rm(list = ls(all.names = TRUE))


# SET RANDOM SEED FOR REPRODUCIBILITY -------------------------------------

set.seed(897231876)  ## 0 - 10,000


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


setDTthreads(threads = 0)  # 0 = all available


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


# STATISTICAL ANALYSIS - LAH - HERITABILITY ---------------------------------------------------

## Run loop to calculate genetic variances for each temperature and population
start <- Sys.time()
if(file.exists("data/Heritability_Bootstrap/Heritability_LAH_var_boot.csv") == FALSE) {
  heritability.tl.boot <- do.call(rbind, mclapply(unique(larval.tl$group), mc.cores = detectCores(), function(grp) {
    ## Filter to only a single temperature treatment
    data.group <- larval.tl %>% filter(group == grp) %>% 
      select(family, dam, sire, block, tl.raw = length_mm)
    
    ## Create a bootstrapped data set from each temperature and group
    bootstrap.data <- do.call(cbind, lapply(1:10000, function(length) {
      ## create a vector of randomly generated data (randomly sample each variable and repeat by nrow(data.temp.group))
      data.family <- do.call(rbind, lapply(unique(data.group$family), function(fam) {
        data.family <- data.group %>% filter(family == fam)
        tl.boot <- sample(data.family$tl.raw, replace = T, size = nrow(data.family))
        data.family.boot <- data.frame(data.family, tl.boot) %>% select(-tl.raw)
      })) %>% 
        mutate(family = factor(family),
               dam = factor(dam),
               sire = factor(sire),
               block = factor(block))
      
      colnames(data.family) <- c(paste0("family", length), paste0("dam", length), paste0("sire", length), 
                                 paste0("block", length), paste0("tl", length))
      data.family
    })) %>% 
      transmute(group = grp, !!!.)
    
    ## Save bootstrapped fish data
    write.csv(bootstrap.data, paste0("data/Heritability_Bootstrap/LAH/Heritability_LAH_Boot_Fish_", grp ,".csv"), row.names = FALSE)
    
    ## Calculate variance components from bootstrapped sample
    bootstrap.data.glmer2 <- resampLmer2(resamp = bootstrap.data, dam = "dam", sire = "sire", response = "tl",
                                          block = "block", start = 1, end = 10000)
    
    bootstrap.data.glmer2.h2 <- data.frame(bootstrap.data.glmer2) %>% 
      mutate(group = grp, rep = 1:n(), trait = "tl") %>% 
      select(group, trait, rep, sire, dam, dam.sire, block, residual = Residual, total = Total, additive, nonadd, maternal)
  }))
  
  ## Save variances for future use
  write.csv(heritability.tl.boot, "data/Heritability_Bootstrap/Heritability_LAH_var_boot.csv", row.names = FALSE)
} else {
  heritability.tl.boot <- fread("data/Heritability_Bootstrap/Heritability_LAH_var_boot.csv")
}

end <- Sys.time()
end - start


# STATISTICAL ANALYSIS - YSV - HERITABILITY ---------------------------------------------------

## Run loop to calculate genetic variances for each temperature and population
start <- Sys.time()
if(file.exists("data/Heritability_Bootstrap/Heritability_YSV_var_boot.csv") == FALSE) {

  heritability.yolk.boot <- do.call(rbind, mclapply(unique(larval.yolk$group), mc.cores = detectCores(), function(grp) {
    ## Filter to only a single temperature treatment
    data.group <- larval.yolk %>% filter(group == grp) %>% 
      select(family, dam, sire, block, yolk.raw = y_vol_mm3)
    
    ## Create a bootstrapped data set from each temperature and group
    bootstrap.data <- do.call(cbind, lapply(1:10000, function(length) {
      ## create a vector of randomly generated data (randomly sample each variable and repeat by nrow(data.temp.group))
      data.family <- do.call(rbind, lapply(unique(data.group$family), function(fam) {
        data.family <- data.group %>% filter(family == fam)
        yolk.boot <- sample(data.family$yolk.raw, replace = T, size = nrow(data.family))
        data.family.boot <- data.frame(data.family, yolk.boot) %>% select(-yolk.raw)
      })) %>% 
        mutate(family = factor(family),
               dam = factor(dam),
               sire = factor(sire),
               block = factor(block))
      
      colnames(data.family) <- c(paste0("family", length), paste0("dam", length), paste0("sire", length), 
                                 paste0("block", length), paste0("yolk", length))
      data.family
    })) %>% 
      transmute(group = grp, !!!.)
    
    ## Save bootstrapped fish data
    write.csv(bootstrap.data, paste0("data/Heritability_Bootstrap/YSV/Heritability_YSV_Boot_Fish_", grp ,".csv"), row.names = FALSE)
    
    bootstrap.data <- fread(grp)
    
    ## Calculate variance components from bootstrapped sample
    bootstrap.data.glmer2 <- resampLmer2(resamp = bootstrap.data, dam = "dam", sire = "sire", response = "yolk",
                                         block = "block", start = 1, end = 10000)
    
    bootstrap.data.glmer2.h2 <- data.frame(bootstrap.data.glmer2) %>% 
      mutate(group = grp, rep = 1:n(), trait = "yolk") %>% 
      select(group, trait, rep, sire, dam, dam.sire, block, residual = Residual, total = Total, additive, nonadd, maternal)
  }))
  
  ## Save variances for future use
  write.csv(heritability.yolk.boot2, "data/Heritability_Bootstrap/Heritability_YSV_Var_Boot.csv", row.names = FALSE)
} else {
  heritability.yolk.boot <- fread("data/Heritability_Bootstrap/Heritability_YSV_Var_Boot.csv")
}

end <- Sys.time()
end - start


# STATISTICAL ANALYSIS - GENERATE OBSERVED HERITABILITY ---------------------------------------

## Length-at-Hatch
heritability.tl.obs <- do.call(rbind, lapply(unique(larval.tl$group), function(grp) {
  ## Filter to only a single temperature treatment
  data.grp <- larval.tl %>% filter(group == grp) %>% 
      select(family, dam, sire, block, length_mm)
    
    obs.tl <- observLmer2(observ = data.grp, dam = "dam", sire = "sire", response = "length_mm", block = "block")
    
    obs.tl.df <- data.frame(group = grp,
                            block = obs.tl$random[4,2],
                            residual = obs.tl$other[1,2],
                            additive = obs.tl$calculation[1,2],
                            nonadd = obs.tl$calculation[2,2],
                            maternal = obs.tl$calculation[3,2]) %>% 
      mutate(pheno = additive + maternal + residual,
             h2.obs = additive / pheno,
             maternal.obs = maternal / pheno) %>% 
      select(group, h2.obs, maternal.obs)
}))

## Yolk-sac Volume
heritability.yolk.obs <- do.call(rbind, lapply(unique(larval.yolk$group), function(grp) {
  ## Filter to only a single temperature treatment
  data.grp <- larval.yolk %>% filter(group == grp) %>% 
    select(family, dam, sire, block, y_vol_mm3)
  
  obs.yolk <- observLmer2(observ = data.grp, dam = "dam", sire = "sire", response = "y_vol_mm3", block = "block")
  
  obs.yolk.df <- data.frame(group = grp,
                            block = obs.yolk$random[4,2],
                            residual = obs.yolk$other[1,2],
                            additive = obs.yolk$calculation[1,2],
                            nonadd = obs.yolk$calculation[2,2],
                            maternal = obs.yolk$calculation[3,2]) %>% 
    mutate(pheno = additive + maternal + residual,
           h2.obs = additive / pheno,
           maternal.obs = maternal / pheno) %>% 
    select(group, h2.obs, maternal.obs)
}))


# CALCULATE THE BIAS-CORRECTED MEAN AND SE FROM BOOTSTRAPPED DISTRIBUTIONS  -------------------

## Length-at-Hatch
heritability.tl.summary <- heritability.tl.boot %>% 
  mutate(pheno = additive + maternal + residual,
         h2.boot = additive / pheno,
         maternal.boot = maternal / pheno) %>% 
  left_join(heritability.tl.obs) %>% 
  group_by(group) %>% 
  summarize(h2.obs = median(h2.obs),
            h2.boot.mean = mean(h2.boot),
            h2.obs.bias = h2.obs - (h2.boot.mean-h2.obs),
            h2.se = sd(h2.boot),
            maternal.obs = median(maternal.obs),
            maternal.boot.mean = mean(maternal.boot),
            maternal.obs.bias = maternal.obs - (maternal.boot.mean-maternal.obs),
            maternal.se = sd(maternal.boot)) %>% 
  mutate(trait = "tl",
         h2.obs.bias = ifelse(h2.obs.bias < 0, 0, h2.obs.bias)) %>% 
  ## Round numeric columns
  mutate_if(is.numeric, round, 2) %>% 
  mutate(population = gsub("_", "-", substr(group, 1, nchar(group)-5)),
         temperature = substr(group, nchar(group)-3, nchar(group)),
         temperature = gsub("_", ".", temperature),
         temperature = as.numeric(gsub("C", "", temperature))) %>% 
  select(group, population, temperature, everything())

## Yolk-sac Volume
heritability.yolk.summary <- heritability.yolk.boot %>% 
  mutate(pheno = additive + maternal + residual,
         h2.boot = additive / pheno,
         maternal.boot = maternal / pheno) %>% 
  left_join(heritability.yolk.obs) %>% 
  group_by(group) %>% 
  summarize(h2.obs = median(h2.obs),
            h2.boot.mean = mean(h2.boot),
            h2.obs.bias = h2.obs - (h2.boot.mean-h2.obs),
            h2.se = sd(h2.boot),
            maternal.obs = median(maternal.obs),
            maternal.boot.mean = mean(maternal.boot),
            maternal.obs.bias = maternal.obs - (maternal.boot.mean-maternal.obs),
            maternal.se = sd(maternal.boot)) %>% 
  mutate(trait = "yolk",
         h2.obs.bias = ifelse(h2.obs.bias < 0, 0, h2.obs.bias)) %>% 
  ## Round numeric columns
  mutate_if(is.numeric, round, 2) %>% 
  mutate(population = gsub("_", "-", substr(group, 1, nchar(group)-5)),
         temperature = substr(group, nchar(group)-3, nchar(group)),
         temperature = gsub("_", ".", temperature),
         temperature = as.numeric(gsub("C", "", temperature))) %>% 
  select(group, population, temperature, everything())


#### CREATE DATA FRAME WITH TEMPERATURE TREATMENTS -----------------------------------------------

temp <- data.frame(group = c("LK_Whitefish_8_0C", "LK_Whitefish_6_9C", "LK_Whitefish_4_0C", "LK_Whitefish_2_2C",
                             "LK_Vendace_8_0C", "LK_Vendace_6_9C", "LK_Vendace_4_0C", "LK_Vendace_2_2C",
                             "LS_Cisco_8_9C", "LS_Cisco_6_9C", "LS_Cisco_4_4C", "LS_Cisco_2_0C",
                             "LO_Cisco_8_9C", "LO_Cisco_6_9C", "LO_Cisco_4_4C", "LO_Cisco_2_0C"),
                   temp.treatment = rep(c("9.0°C", "7.0°C", "4.5°C", "2.0°C"), 4))


# COMBINE ALL TRAITS --------------------------------------------------------------------------

heritability.all <- bind_rows(heritability.tl.summary, heritability.yolk.summary) %>%
  left_join(temp) %>% 
  mutate(trait = factor(trait, ordered = TRUE, levels = c("tl", "yolk"),
                        labels = c("LAH", "YSV")),
         temp.treatment = factor(temp.treatment, ordered = TRUE, levels = c("2.0°C", "4.5°C", "7.0°C", "9.0°C")),
         population = factor(population, ordered = TRUE, levels = c("LK-Vendace", "LK-Whitefish", "LS-Cisco", "LO-Cisco")),
         h2.obs.bias = ifelse(h2.obs.bias > 1, 1, h2.obs.bias),
         maternal.obs.bias = ifelse(maternal.obs.bias < 0, 0, maternal.obs.bias)) %>% 
  filter(population != "LK-Whitefish")


#### VISUALIZATION - HERITABILITY --------------------------------------------

## Heritability
plot.h2.lah <- ggplot(filter(heritability.all, trait == "LAH"), aes(x = temperature, y = (h2.obs.bias * 100), group = population, color = population, shape = population, linetype = population)) + 
  geom_line(size = 1.0, position = position_dodge(0.13)) +
  geom_point(size = 5, position = position_dodge(0.13)) +
  annotate("text", label = "A", x = 2.0, y = 65, size = 7) +
  geom_errorbar(aes(ymin = ifelse((h2.obs.bias - h2.se) * 100 < 0, 0, (h2.obs.bias - h2.se) * 100), 
                    ymax = ifelse((h2.obs.bias + h2.se) * 100 > 100, 100, (h2.obs.bias + h2.se) * 100)), 
                position = position_dodge(0.13),
                size = 1.0, width = 0.25, linetype = "solid", show.legend = FALSE) +
  #scale_color_manual("combine", values = c("#000000", "#717171" ,"#9f9e9f", "#c6c5c6"),
  #labels = c("LK-Vendace   ", "LK-Whitefish   ", "LS-Cisco   ", "LO-Cisco")) +
  scale_color_manual("combine", values = c("#000000","#9f9e9f", "#c6c5c6"),
                     labels = c("LK-Vendace   ", "LS-Cisco   ", "LO-Cisco")) +
  #scale_shape_manual("combine", values = c(2, 5, 1, 0), 
  #labels = c("LK-Vendace   ", "LK-Whitefish   ", "LS-Cisco   ", "LO-Cisco")) +
  scale_shape_manual("combine", values = c(2, 1, 0), 
                     labels = c("LK-Vendace   ", "LS-Cisco   ", "LO-Cisco")) +
  #scale_linetype_manual("combine", values = c("solid", "dashed", "dotted", "solid"), 
  #labels = c("LK-Vendace   ", "LK-Whitefish   ", "LS-Cisco   ", "LO-Cisco")) +
  scale_linetype_manual("combine", values = c("solid", "dotted", "solid"), 
                        labels = c("LK-Vendace   ", "LS-Cisco   ", "LO-Cisco")) +
  scale_y_continuous(limits = c(-2, 70), breaks = seq(0, 100, 10), expand = c(0, 0)) +
  scale_x_continuous(limits = c(1.75, 9.15), breaks = c(2, 4, 4.4, 6.9, 8, 8.9), expand = c(0, 0)) +
  labs(x = "Incubation Temperature (°C)", y = "Narrow-sense Heritability (%)") +
  theme_bw() +
  theme(axis.title.x = element_text(color = "Black", size = 22, margin = margin(15, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 22, margin = margin(0, 15, 0, 0)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 20),
        legend.key.size = unit(1.25, 'cm'),
        legend.position = "top",
        plot.margin = unit(c(5, 5, 5, 5), 'mm'))


plot.h2.ysv <- ggplot(filter(heritability.all, trait == "YSV"), aes(x = temperature, y = (h2.obs.bias * 100), group = population, color = population, shape = population, linetype = population)) + 
  geom_line(size = 1.0, position = position_dodge(0.13)) +
  geom_point(size = 5, position = position_dodge(0.13)) +
  annotate("text", label = "B", x = 2.0, y = 46, size = 7) +
  geom_errorbar(aes(ymin = ifelse((h2.obs.bias - h2.se) * 100 < 0, 0, (h2.obs.bias - h2.se) * 100), 
                    ymax = ifelse((h2.obs.bias + h2.se) * 100 > 100, 100, (h2.obs.bias + h2.se) * 100)), 
                position = position_dodge(0.13),
                size = 1.0, width = 0.25, linetype = "solid", show.legend = FALSE) +
  #scale_color_manual("combine", values = c("#000000", "#717171" ,"#9f9e9f", "#c6c5c6"),
  #labels = c("LK-Vendace   ", "LK-Whitefish   ", "LS-Cisco   ", "LO-Cisco")) +
  scale_color_manual("combine", values = c("#000000","#9f9e9f", "#c6c5c6"),
                     labels = c("LK-Vendace   ", "LS-Cisco   ", "LO-Cisco")) +
  #scale_shape_manual("combine", values = c(2, 5, 1, 0), 
  #labels = c("LK-Vendace   ", "LK-Whitefish   ", "LS-Cisco   ", "LO-Cisco")) +
  scale_shape_manual("combine", values = c(2, 1, 0), 
                     labels = c("LK-Vendace   ", "LS-Cisco   ", "LO-Cisco")) +
  #scale_linetype_manual("combine", values = c("solid", "dashed", "dotted", "solid"), 
  #labels = c("LK-Vendace   ", "LK-Whitefish   ", "LS-Cisco   ", "LO-Cisco")) +
  scale_linetype_manual("combine", values = c("solid", "dotted", "solid"), 
                        labels = c("LK-Vendace   ", "LS-Cisco   ", "LO-Cisco")) +
  scale_y_continuous(limits = c(-2, 50), breaks = seq(0, 100, 10), expand = c(0, 0)) +
  scale_x_continuous(limits = c(1.75, 9.15), breaks = c(2, 4, 4.4, 6.9, 8, 8.9), expand = c(0, 0)) +
  labs(x = "Incubation Temperature (°C)", y = "Narrow-sense Heritability (%)") +
  theme_bw() +
  theme(axis.title.x = element_text(color = "Black", size = 22, margin = margin(15, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 22, margin = margin(0, 15, 0, 0)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 20),
        legend.key.size = unit(1.25, 'cm'),
        legend.position = "top",
        plot.margin = unit(c(5, 5, 5, 5), 'mm'))

## Combine all figures
plot.h2.all <- grid.arrange(arrangeGrob(textGrob(""), 
                                     get_legend(plot.h2.lah),
                                     nrow = 1,
                                     widths = c(0.09, 1)),
                         arrangeGrob(plot.h2.lah + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank()),
                                     plot.h2.ysv + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank()),
                                     nrow = 2,
                                     left = textGrob("Narrow-sense Heritability (%)", y = 0.52, rot = 90, gp = gpar(cex = 1.75, fontfamily = "Arial")),
                                     bottom = textGrob("Mean Incubation Temperature (°C)", x = 0.545, gp = gpar(cex = 1.75, fontfamily = "Arial"))),
                         heights = c(0.025, 1)
)

ggsave("figures/2020-Larval-Heritability-SE-Line.png", plot = plot.h2.all, width = 11, height = 10, dpi = 300)

 ## Heritability
plot.m2.lah <- ggplot(filter(heritability.all, trait == "LAH"), aes(x = temperature, y = (maternal.obs.bias * 100), group = population, color = population, shape = population, linetype = population)) + 
  geom_line(size = 1.0, position = position_dodge(0.13)) +
  geom_point(size = 5, position = position_dodge(0.13)) +
  annotate("text", label = "A", x = 2.0, y = 95, size = 7) +
  geom_errorbar(aes(ymin = ifelse((maternal.obs.bias - maternal.se) * 100 < 0, 0, (maternal.obs.bias - maternal.se) * 100), 
                    ymax = ifelse((maternal.obs.bias + maternal.se) * 100 > 100, 100, (maternal.obs.bias + maternal.se) * 100)), 
                position = position_dodge(0.13),
                size = 1.0, width = 0.25, linetype = "solid", show.legend = FALSE) +
  #scale_color_manual("combine", values = c("#000000", "#717171" ,"#9f9e9f", "#c6c5c6"),
  #labels = c("LK-Vendace   ", "LK-Whitefish   ", "LS-Cisco   ", "LO-Cisco")) +
  scale_color_manual("combine", values = c("#000000","#9f9e9f", "#c6c5c6"),
                     labels = c("LK-Vendace   ", "LS-Cisco   ", "LO-Cisco")) +
  #scale_shape_manual("combine", values = c(2, 5, 1, 0), 
  #labels = c("LK-Vendace   ", "LK-Whitefish   ", "LS-Cisco   ", "LO-Cisco")) +
  scale_shape_manual("combine", values = c(2, 1, 0), 
                     labels = c("LK-Vendace   ", "LS-Cisco   ", "LO-Cisco")) +
  #scale_linetype_manual("combine", values = c("solid", "dashed", "dotted", "solid"), 
  #labels = c("LK-Vendace   ", "LK-Whitefish   ", "LS-Cisco   ", "LO-Cisco")) +
  scale_linetype_manual("combine", values = c("solid", "dotted", "solid"), 
                        labels = c("LK-Vendace   ", "LS-Cisco   ", "LO-Cisco")) +
  scale_y_continuous(limits = c(-2, 102), breaks = seq(0, 100, 20), expand = c(0, 0)) +
  scale_x_continuous(limits = c(1.75, 9.15), breaks = c(2, 4, 4.4, 6.9, 8, 8.9), expand = c(0, 0)) +
  labs(x = "Incubation Temperature (°C)", y = "Maternal Effect (%)") +
  theme_bw() +
  theme(axis.title.x = element_text(color = "Black", size = 22, margin = margin(15, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 22, margin = margin(0, 15, 0, 0)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 20),
        legend.key.size = unit(1.25, 'cm'),
        legend.position = "top",
        plot.margin = unit(c(5, 5, 5, 5), 'mm'))


plot.m2.ysv <- ggplot(filter(heritability.all, trait == "YSV"), aes(x = temperature, y = (maternal.obs.bias * 100), group = population, color = population, shape = population, linetype = population)) + 
  geom_line(size = 1.0, position = position_dodge(0.13)) +
  geom_point(size = 5, position = position_dodge(0.13)) +
  annotate("text", label = "B", x = 2.0, y = 65, size = 7) +
  geom_errorbar(aes(ymin = ifelse((maternal.obs.bias - maternal.se) * 100 < 0, 0, (maternal.obs.bias - maternal.se) * 100), 
                    ymax = ifelse((maternal.obs.bias + maternal.se) * 100 > 100, 100, (maternal.obs.bias + maternal.se) * 100)), 
                position = position_dodge(0.13),
                size = 1.0, width = 0.25, linetype = "solid", show.legend = FALSE) +
  #scale_color_manual("combine", values = c("#000000", "#717171" ,"#9f9e9f", "#c6c5c6"),
  #labels = c("LK-Vendace   ", "LK-Whitefish   ", "LS-Cisco   ", "LO-Cisco")) +
  scale_color_manual("combine", values = c("#000000","#9f9e9f", "#c6c5c6"),
                     labels = c("LK-Vendace   ", "LS-Cisco   ", "LO-Cisco")) +
  #scale_shape_manual("combine", values = c(2, 5, 1, 0), 
  #labels = c("LK-Vendace   ", "LK-Whitefish   ", "LS-Cisco   ", "LO-Cisco")) +
  scale_shape_manual("combine", values = c(2, 1, 0), 
                     labels = c("LK-Vendace   ", "LS-Cisco   ", "LO-Cisco")) +
  #scale_linetype_manual("combine", values = c("solid", "dashed", "dotted", "solid"), 
  #labels = c("LK-Vendace   ", "LK-Whitefish   ", "LS-Cisco   ", "LO-Cisco")) +
  scale_linetype_manual("combine", values = c("solid", "dotted", "solid"), 
                        labels = c("LK-Vendace   ", "LS-Cisco   ", "LO-Cisco")) +
  scale_y_continuous(limits = c(-2, 70), breaks = seq(0, 100, 10), expand = c(0, 0)) +
  scale_x_continuous(limits = c(1.75, 9.15), breaks = c(2, 4, 4.4, 6.9, 8, 8.9), expand = c(0, 0)) +
  labs(x = "Incubation Temperature (°C)", y = "Maternal Effect (%)") +
  theme_bw() +
  theme(axis.title.x = element_text(color = "Black", size = 22, margin = margin(15, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 22, margin = margin(0, 15, 0, 0)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 20),
        legend.key.size = unit(1.25, 'cm'),
        legend.position = "top",
        plot.margin = unit(c(5, 5, 5, 5), 'mm'))

## Combine all figures
plot.m2.all <- grid.arrange(arrangeGrob(textGrob(""), 
                                        get_legend(plot.m2.lah),
                                        nrow = 1,
                                        widths = c(0.09, 1)),
                            arrangeGrob(plot.m2.lah + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank()),
                                        plot.m2.ysv + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank()),
                                        nrow = 2,
                                        left = textGrob("Maternal Effect (%)", y = 0.52, rot = 90, gp = gpar(cex = 1.75, fontfamily = "Arial")),
                                        bottom = textGrob("Mean Incubation Temperature (°C)", x = 0.545, gp = gpar(cex = 1.75, fontfamily = "Arial"))),
                            heights = c(0.025, 1)
)

ggsave("figures/2020-Larval-Maternal-SE-Line.png", plot = plot.m2.all, width = 11, height = 10, dpi = 300)
