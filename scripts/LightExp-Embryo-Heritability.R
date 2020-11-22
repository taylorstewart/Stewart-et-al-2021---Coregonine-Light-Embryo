#### CLEAR THE ENVIRONMENT FIRST -----------------------------------------------------------------

rm(list = ls(all.names = TRUE))


#### SET RANDOM SEED FOR REPRODUCIBILITY ---------------------------------------------------------

set.seed(897231876)  ## 0 - 10,000


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

setDTthreads(threads = 0)  # 0 = all available


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


#### STATISTICAL ANALYSIS - SURVIVAL - HERITABILITY ----------------------------------------------

## Run loop to calculate genetic variances for each temperature and population
start <- Sys.time()
if(file.exists("data/Heritability_Bootstrap/Heritability_ES_Var_Boot.csv") == FALSE) {
  heritability.survival.boot <- do.call(rbind, mclapply(unique(hatch.survival$group), mc.cores = detectCores(), function(grp) {
    ## Filter to only a single temperature treatment
    data.group <- hatch.survival %>% filter(group == grp) %>% 
      select(family, dam, sire, block, hatch.raw = hatch)
    
    ## Create a bootstrapped data set from each temperature and group
    bootstrap.data <- do.call(cbind, lapply(1:10000, function(length) {
      ## create a vector of randomly generated data (randomly sample each variable and repeat by nrow(data.group))
      repeat {
        data.family <- do.call(rbind, lapply(unique(data.group$family), function(fam) {
          data.family <- data.group %>% filter(family == fam)
          hatch.boot <- sample(data.family$hatch.raw, replace = T, size = nrow(data.family))
          data.family.boot <- data.frame(data.family, hatch.boot) %>% select(-hatch.raw)
        })) %>% 
          mutate(family = factor(family),
                 dam = factor(dam),
                 sire = factor(sire),
                 block = factor(block))
        if(sum(data.family$hatch.boot) != 1402) {
          break
        }
        }
      
      colnames(data.family) <- c(paste0("family", length), paste0("dam", length), paste0("sire", length), 
                                 paste0("block", length), paste0("hatch", length))
      data.family
    })) %>% 
      transmute(group = grp, !!!.)
    
    ## Save bootstrapped fish data
    write.csv(bootstrap.data, paste0("data/Heritability_Bootstrap/ES/Heritability_ES_Boot_Fish_", grp, ".csv"), row.names = FALSE)
    
    ## Calculate variance components from bootstrapped sample
    bootstrap.data.glmer2 <- resampGlmer2(resamp = bootstrap.data, dam = "dam", sire = "sire", response = "hatch",
                                          block = "block", fam_link = binomial(logit), start = 1, end = 10000)
    
    bootstrap.data.glmer2.h2 <- data.frame(bootstrap.data.glmer2) %>% 
      mutate(group = grp, rep = 1:n(), trait = "survival") %>% 
      select(group, trait, rep, sire, dam, dam.sire, block, residual = Residual, total = Total, additive, nonadd, maternal)
  }))
  
  ## Save variances for future use
  write.csv(heritability.survival.boot, "data/Heritability_Bootstrap/Heritability_ES_Var_Boot.csv", row.names = FALSE)
} else {
  heritability.survival.boot <- fread("data/Heritability_Bootstrap/Heritability_ES_Var_Boot.csv")
}

end <- Sys.time()
end - start


#### STATISTICAL ANALYSIS - INCUBATION PERIOD (DPF) - HERITABILITY -------------------------------

## Run loop to calculate genetic variances for each temperature and population
start <- Sys.time()
if(file.exists("data/Heritability_Bootstrap/Heritability_DPF_Var_Boot.csv") == FALSE) {
  heritability.dpf.boot <- do.call(rbind, mclapply(unique(hatch.dpf$group), mc.cores = detectCores(), function(grp) {
    ## Filter to only a single temperature treatment
    data.group <- hatch.dpf %>% filter(group == grp) %>% 
      select(family, dam, sire, block, dpf.raw = dpf)
    
    ## Create a bootstrapped data set from each temperature and group
    bootstrap.data <- do.call(cbind, lapply(1:10000, function(length) {
      ## create a vector of randomly generated data (randomly sample each variable and repeat by nrow(data.group))
      data.family <- do.call(rbind, lapply(unique(data.group$family), function(fam) {
        data.family <- data.group %>% filter(family == fam)
        dpf.boot <- sample(data.family$dpf.raw, replace = T, size = nrow(data.family))
        data.family.boot <- data.frame(data.family, dpf.boot) %>% select(-dpf.raw)
      })) %>% 
        mutate(family = factor(family),
               dam = factor(dam),
               sire = factor(sire),
               block = factor(block))
      
      colnames(data.family) <- c(paste0("family", length), paste0("dam", length), paste0("sire", length), 
                                 paste0("block", length), paste0("dpf", length))
      data.family
    })) %>% 
      transmute(group = grp, !!!.)
    
    ## Save bootstrapped fish data
    write.csv(bootstrap.data, paste0("data/Heritability_Bootstrap/DPF/Heritability_DPF_Boot_Fish_", grp ,".csv"), row.names = FALSE)
    
    ## Calculate variance components from bootstrapped sample
    bootstrap.data.glmer2 <- resampLmer2(resamp = bootstrap.data, dam = "dam", sire = "sire", response = "dpf",
                                          block = "block", start = 1, end = 10000)
    
    bootstrap.data.glmer2.h2 <- data.frame(bootstrap.data.glmer2) %>% 
      mutate(group = grp, rep = 1:n(), trait = "dpf") %>% 
      select(group, trait, rep, sire, dam, dam.sire, block, residual = Residual, total = Total, additive, nonadd, maternal)
  }))
  
  ## Save variances for future use
  write.csv(heritability.dpf.boot, "data/Heritability_Bootstrap/Heritability_DPF_Var_Boot.csv", row.names = FALSE)
} else {
  heritability.dpf.boot <- fread("data/Heritability_Bootstrap/Heritability_DPF_Var_Boot.csv")
}

end <- Sys.time()
end - start


#### STATISTICAL ANALYSIS - INCUBATION PERIOD (ADD) - HERITABILITY -------------------------------

## Run loop to calculate genetic variances for each temperature and population
start <- Sys.time()
if(file.exists("data/Heritability_Bootstrap/Heritability_ADD_Var_Boot.csv") == FALSE) {
  heritability.ADD.boot <- do.call(rbind, mclapply(unique(hatch.ADD$group), mc.cores = detectCores(), function(grp) {
    ## Filter to only a single temperature treatment
    data.group <- hatch.ADD %>% filter(group == grp) %>% 
      select(family, dam, sire, block, ADD.raw = ADD)
    
    ## Create a bootstrapped data set from each temperature and group
    bootstrap.data <- do.call(cbind, lapply(1:10000, function(length) {
      ## create a vector of randomly generated data (randomly sample each variable and repeat by nrow(data.group))
      data.family <- do.call(rbind, lapply(unique(data.group$family), function(fam) {
        data.family <- data.group %>% filter(family == fam)
        ADD.boot <- sample(data.family$ADD.raw, replace = T, size = nrow(data.family))
        data.family.boot <- data.frame(data.family, ADD.boot) %>% select(-ADD.raw)
      })) %>% 
        mutate(family = factor(family),
               dam = factor(dam),
               sire = factor(sire),
               block = factor(block))
      
      colnames(data.family) <- c(paste0("family", length), paste0("dam", length), paste0("sire", length), 
                                 paste0("block", length), paste0("ADD", length))
      data.family
    })) %>% 
      transmute(group = grp, !!!.)
    
    ## Save bootstrapped fish data
    write.csv(bootstrap.data, paste0("data/Heritability_Bootstrap/ADD/Heritability_ADD_Boot_Fish_", grp ,".csv"), row.names = FALSE)
    
    ## Calculate variance components from bootstrapped sample
    bootstrap.data.glmer2 <- resampLmer2(resamp = bootstrap.data, dam = "dam", sire = "sire", response = "ADD",
                                          block = "block", start = 1, end = 10000)
    
    bootstrap.data.glmer2.h2 <- data.frame(bootstrap.data.glmer2) %>% 
      mutate(group = grp, rep = 1:n(), trait = "ADD") %>% 
      select(group, trait, rep, sire, dam, dam.sire, block, residual = Residual, total = Total, additive, nonadd, maternal)
  }))
  
  ## Save variances for future use
  write.csv(heritability.ADD.boot, "data/Heritability_Bootstrap/Heritability_ADD_Var_Boot.csv", row.names = FALSE)
} else {
  heritability.ADD.boot <- fread("data/Heritability_Bootstrap/Heritability_ADD_Var_Boot.csv")
}

end <- Sys.time()
end - start


#### STATISTICAL ANALYSIS - GENERATE OBSERVED HERITABILITY ---------------------------------------

## Embryo Survival
heritability.survival.obs <- do.call(rbind, lapply(unique(hatch.survival$group), function(grp) {
  ## Filter to only a single group
  data.group <- hatch.survival %>% filter(group == grp) %>% 
      select(family, dam, sire, block, hatch)
  
  obs.survival <- observGlmer2(observ = data.group, dam = "dam", sire = "sire", response = "hatch",
                               block = "block", fam_link = binomial(logit))
  
  obs.survival.df <- data.frame(group = grp,
                                block = obs.survival$random[4,2],
                                residual = obs.survival$other[1,2],
                                additive = obs.survival$calculation[1,2],
                                nonadd = obs.survival$calculation[2,2],
                                maternal = obs.survival$calculation[3,2]) %>% 
    mutate(pheno = additive + maternal + residual,
           h2.obs = additive / pheno,
           maternal.obs = maternal / pheno) %>% 
    select(group, h2.obs, maternal.obs)
}))

## DPF
heritability.dpf.obs <- do.call(rbind, lapply(unique(hatch.dpf$group), function(grp) {
  ## Filter to only a single group
  data.group <- hatch.dpf %>% filter(group == grp) %>% 
      select(family, dam, sire, block, dpf)
    
  obs.dpf <- observLmer2(observ = data.group, dam = "dam", sire = "sire", response = "dpf", block = "block")
  
  obs.dpf.df <- data.frame(group = grp,
                           block = obs.dpf$random[4,2],
                           residual = obs.dpf$other[1,2],
                           additive = obs.dpf$calculation[1,2],
                           nonadd = obs.dpf$calculation[2,2],
                           maternal = obs.dpf$calculation[3,2]) %>% 
    mutate(pheno = additive + maternal + residual,
           h2.obs = additive / pheno,
           maternal.obs = maternal / pheno) %>% 
    select(group, h2.obs, maternal.obs)
}))

## ADD
heritability.ADD.obs <- do.call(rbind, lapply(unique(hatch.ADD$group), function(grp) {
  ## Filter to only a single group
  data.group <- hatch.ADD %>% filter(group == grp) %>% 
      select(family, dam, sire, block, ADD)
  
  obs.ADD <- observLmer2(observ = data.group, dam = "dam", sire = "sire", response = "ADD", block = "block")
  
  obs.ADD.df <- data.frame(group = grp,
                           block = obs.ADD$random[4,2],
                           residual = obs.ADD$other[1,2],
                           additive = obs.ADD$calculation[1,2],
                           nonadd = obs.ADD$calculation[2,2],
                           maternal = obs.ADD$calculation[3,2]) %>% 
    mutate(pheno = additive + maternal + residual,
           h2.obs = additive / pheno,
           maternal.obs = maternal / pheno) %>% 
    select(group, h2.obs, maternal.obs)
}))


#### CALCULATE THE BIAS-CORRECTED MEAN AND SE FROM BOOTSTRAPPED DISTRIBUTIONS  -------------------

## Embryo Survival
heritability.survival.summary <- heritability.survival.boot %>% 
  mutate(pheno = additive + maternal + residual,
         h2.boot = additive / pheno,
         maternal.boot = maternal / pheno) %>% 
  left_join(heritability.survival.obs) %>% 
  group_by(group) %>% 
  summarize(h2.obs = median(h2.obs),
            h2.boot.mean = mean(h2.boot),
            h2.obs.bias = h2.obs - (h2.boot.mean-h2.obs),
            h2.se = sd(h2.boot),
            maternal.obs = median(maternal.obs),
            maternal.boot.mean = mean(maternal.boot),
            maternal.obs.bias = maternal.obs - (maternal.boot.mean-maternal.obs),
            maternal.se = sd(maternal.boot)) %>% 
  mutate(trait = "survival",
         h2.obs.bias = ifelse(h2.obs.bias < 0, 0, h2.obs.bias)) %>% 
  ## Round numeric columns
  mutate_if(is.numeric, round, 2) %>% 
  mutate(population = gsub("_", "", substr(group, 1, 8)),
         treatment = gsub("_", "", substr(group, 9, nchar(group)))) %>% 
  select(group, population, treatment, everything())

## DPF
heritability.dpf.summary <- heritability.dpf.boot %>% 
  mutate(pheno = additive + maternal + residual,
         h2.boot = additive / pheno,
         maternal.boot = maternal / pheno) %>% 
  left_join(heritability.dpf.obs) %>% 
  group_by(group) %>% 
  summarize(h2.obs = median(h2.obs),
            h2.boot.mean = mean(h2.boot),
            h2.obs.bias = h2.obs - (h2.boot.mean-h2.obs),
            h2.se = sd(h2.boot),
            maternal.obs = median(maternal.obs),
            maternal.boot.mean = mean(maternal.boot),
            maternal.obs.bias = maternal.obs - (maternal.boot.mean-maternal.obs),
            maternal.se = sd(maternal.boot)) %>% 
  mutate(trait = "dpf",
         h2.obs.bias = ifelse(h2.obs.bias < 0, 0, h2.obs.bias)) %>% 
  ## Round numeric columns
  mutate_if(is.numeric, round, 2) %>% 
  mutate(population = gsub("_", "", substr(group, 1, 8)),
         treatment = gsub("_", "", substr(group, 9, nchar(group)))) %>% 
  select(group, population, treatment, everything())

## ADD
heritability.ADD.summary <- heritability.ADD.boot %>% 
  mutate(pheno = additive + maternal + residual,
         h2.boot = additive / pheno,
         maternal.boot = maternal / pheno) %>% 
  left_join(heritability.ADD.obs) %>% 
  group_by(group) %>% 
  summarize(h2.obs = median(h2.obs),
            h2.boot.mean = mean(h2.boot),
            h2.obs.bias = h2.obs - (h2.boot.mean-h2.obs),
            h2.se = sd(h2.boot),
            maternal.obs = median(maternal.obs),
            maternal.boot.mean = mean(maternal.boot),
            maternal.obs.bias = maternal.obs - (maternal.boot.mean-maternal.obs),
            maternal.se = sd(maternal.boot)) %>% 
  mutate(trait = "ADD",
         h2.obs.bias = ifelse(h2.obs.bias < 0, 0, h2.obs.bias)) %>% 
  ## Round numeric columns
  mutate_if(is.numeric, round, 2) %>% 
  mutate(population = gsub("_", "", substr(group, 1, 8)),
         treatment = gsub("_", "", substr(group, 9, nchar(group)))) %>% 
  select(group, population, treatment, everything())


#### COMBINE ALL TRAITS --------------------------------------------------------------------------

heritability.all <- bind_rows(heritability.survival.summary, heritability.ADD.summary, heritability.dpf.summary) %>% 
  mutate(trait = factor(trait, ordered = TRUE, levels = c("survival", "dpf", "ADD"), 
                        labels = c("ES", "DPF", "ADD")),
         population = factor(population, ordered = TRUE, levels = c("Superior", "Ontario")),
         treatment = factor(treatment, ordered = TRUE, levels = c("High", "Medium", "Low")),
         h2.obs.bias = ifelse(h2.obs.bias > 1, 1, h2.obs.bias),
         maternal.obs.bias = ifelse(maternal.obs.bias < 0, 0, maternal.obs.bias))


#### VISUALIZATION - HERITABILITY ----------------------------------------------------------------

## Heritability
#plot.h2.es <- 
ggplot(filter(heritability.all, trait == "ES"), aes(x = treatment, y = (h2.obs.bias * 100), group = population, color = population, shape = population, linetype = population)) + 
  geom_line(size = 1.0, position = position_dodge(0.13)) +
  geom_point(size = 5, position = position_dodge(0.13)) +
  annotate("text", label = "A", x = 1.0, y = 46, size = 7) +
  geom_errorbar(aes(ymin = ifelse((h2.obs.bias - h2.se) * 100 < 0, 0, (h2.obs.bias - h2.se) * 100), 
                    ymax = ifelse((h2.obs.bias + h2.se) * 100 > 100, 100, (h2.obs.bias + h2.se) * 100)), 
                position = position_dodge(0.13),
                size = 1.0, width = 0.25, linetype = "solid", show.legend = FALSE) +
  scale_color_manual("combine", values = c("#000000", "#717171"),
                     labels = c("Superior   ", "Ontario   ")) +
  scale_shape_manual("combine", values = c(2, 5), 
                     labels = c("Superior   ", "Ontario   ")) +
  scale_linetype_manual("combine", values = c("solid", "dashed"), 
                        labels = c("Superior   ", "Ontario   ")) +
  scale_y_continuous(limits = c(-2, 50), breaks = seq(0, 70, 10), expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0.1)) +
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


#plot.h2.dpf <- 
ggplot(filter(heritability.all, trait == "DPF"), aes(x = treatment, y = (h2.obs.bias * 100), group = population, color = population, shape = population, linetype = population)) + 
  geom_line(size = 1.0, position = position_dodge(0.13)) +
  geom_point(size = 5, position = position_dodge(0.13)) +
  annotate("text", label = "B", x = 1.0, y = 46, size = 7) +
  geom_errorbar(aes(ymin = ifelse((h2.obs.bias - h2.se) * 100 < 0, 0, (h2.obs.bias - h2.se) * 100), 
                    ymax = ifelse((h2.obs.bias + h2.se) * 100 > 100, 100, (h2.obs.bias + h2.se) * 100)), 
                position = position_dodge(0.13),
                size = 1.0, width = 0.25, linetype = "solid", show.legend = FALSE) +
  scale_color_manual("combine", values = c("#000000", "#717171"),
                     labels = c("Superior   ", "Ontario   ")) +
  scale_shape_manual("combine", values = c(2, 5), 
                     labels = c("Superior   ", "Ontario   ")) +
  scale_linetype_manual("combine", values = c("solid", "dashed"), 
                        labels = c("Superior   ", "Ontario   ")) +
  scale_y_continuous(limits = c(-2, 50), breaks = seq(0, 70, 10), expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0.1)) +
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

#plot.h2.add <- 
ggplot(filter(heritability.all, trait == "ADD"), aes(x = treatment, y = (h2.obs.bias * 100), group = population, color = population, shape = population, linetype = population)) + 
  geom_line(size = 1.0, position = position_dodge(0.13)) +
  geom_point(size = 5, position = position_dodge(0.13)) +
  annotate("text", label = "C", x = 1.0, y = 46, size = 7) +
  geom_errorbar(aes(ymin = ifelse((h2.obs.bias - h2.se) * 100 < 0, 0, (h2.obs.bias - h2.se) * 100), 
                    ymax = ifelse((h2.obs.bias + h2.se) * 100 > 100, 100, (h2.obs.bias + h2.se) * 100)), 
                position = position_dodge(0.13),
                size = 1.0, width = 0.25, linetype = "solid", show.legend = FALSE) +
  scale_color_manual("combine", values = c("#000000", "#717171"),
                     labels = c("Superior   ", "Ontario   ")) +
  scale_shape_manual("combine", values = c(2, 5), 
                     labels = c("Superior   ", "Ontario   ")) +
  scale_linetype_manual("combine", values = c("solid", "dashed"), 
                        labels = c("Superior   ", "Ontario   ")) +
  scale_y_continuous(limits = c(-2, 50), breaks = seq(0, 70, 10), expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0.1)) +
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
                                     get_legend(plot.h2.es),
                                     nrow = 1,
                                     widths = c(0.09, 1)),
                         arrangeGrob(plot.h2.es + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank()),
                                     plot.h2.dpf + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank()),
                                     plot.h2.add + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank()),
                                     nrow = 3,
                                     left = textGrob("Narrow-sense Heritability (%)", y = 0.52, rot = 90, gp = gpar(cex = 1.75, fontfamily = "Arial")),
                                     bottom = textGrob("Incubation Temperature Treatment (°C)", x = 0.545, gp = gpar(cex = 1.75, fontfamily = "Arial"))),
                         heights = c(0.025, 1)
)

ggsave("figures/2020-Embryo-Heritability-SE-Line.png", plot = plot.h2.all, width = 11, height = 15, dpi = 300)


## Maternal Effects
#plot.m2.es <- 
ggplot(filter(heritability.all, trait == "ES"), aes(x = treatment, y = (maternal.obs.bias * 100), group = population, color = population, shape = population, linetype = population)) + 
  geom_line(size = 1.0, position = position_dodge(0.13)) +
  geom_point(size = 5, position = position_dodge(0.13)) +
  annotate("text", label = "A", x = 1.0, y = 55, size = 7) +
  geom_errorbar(aes(ymin = ifelse((maternal.obs.bias - maternal.se) * 100 < 0, 0, (maternal.obs.bias - maternal.se) * 100), 
                    ymax = ifelse((maternal.obs.bias + maternal.se) * 100 > 100, 100, (maternal.obs.bias + maternal.se) * 100)), 
                position = position_dodge(0.13),
                size = 1.0, width = 0.25, linetype = "solid", show.legend = FALSE) +
  scale_color_manual("combine", values = c("#000000", "#717171"),
                     labels = c("Superior   ", "Ontario   ")) +
  scale_shape_manual("combine", values = c(2, 5), 
                     labels = c("Superior   ", "Ontario   ")) +
  scale_linetype_manual("combine", values = c("solid", "dashed"), 
                        labels = c("Superior   ", "Ontario   ")) +
  scale_y_continuous(limits = c(-2, 60), breaks = seq(0, 70, 10), expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0.1)) +
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


#plot.m2.dpf <- 
ggplot(filter(heritability.all, trait == "DPF"), aes(x = treatment, y = (maternal.obs.bias * 100), group = population, color = population, shape = population, linetype = population)) + 
  geom_line(size = 1.0, position = position_dodge(0.13)) +
  geom_point(size = 5, position = position_dodge(0.13)) +
  annotate("text", label = "B", x = 1.0, y = 55, size = 7) +
  geom_errorbar(aes(ymin = ifelse((maternal.obs.bias - maternal.se) * 100 < 0, 0, (maternal.obs.bias - maternal.se) * 100), 
                    ymax = ifelse((maternal.obs.bias + maternal.se) * 100 > 100, 100, (maternal.obs.bias + maternal.se) * 100)), 
                position = position_dodge(0.13),
                size = 1.0, width = 0.25, linetype = "solid", show.legend = FALSE) +
  scale_color_manual("combine", values = c("#000000", "#717171"),
                     labels = c("Superior   ", "Ontario   ")) +
  scale_shape_manual("combine", values = c(2, 5), 
                     labels = c("Superior   ", "Ontario   ")) +
  scale_linetype_manual("combine", values = c("solid", "dashed"), 
                        labels = c("Superior   ", "Ontario   ")) +
  scale_y_continuous(limits = c(-2, 60), breaks = seq(0, 70, 10), expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0.1)) +
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

#plot.m2.add <- 
ggplot(filter(heritability.all, trait == "ADD"), aes(x = treatment, y = (maternal.obs.bias * 100), group = population, color = population, shape = population, linetype = population)) + 
  geom_line(size = 1.0, position = position_dodge(0.13)) +
  geom_point(size = 5, position = position_dodge(0.13)) +
  annotate("text", label = "C", x = 1.0, y = 55, size = 7) +
  geom_errorbar(aes(ymin = ifelse((maternal.obs.bias - maternal.se) * 100 < 0, 0, (maternal.obs.bias - maternal.se) * 100), 
                    ymax = ifelse((maternal.obs.bias + maternal.se) * 100 > 100, 100, (maternal.obs.bias + maternal.se) * 100)), 
                position = position_dodge(0.13),
                size = 1.0, width = 0.25, linetype = "solid", show.legend = FALSE) +
  scale_color_manual("combine", values = c("#000000", "#717171"),
                     labels = c("Superior   ", "Ontario   ")) +
  scale_shape_manual("combine", values = c(2, 5), 
                     labels = c("Superior   ", "Ontario   ")) +
  scale_linetype_manual("combine", values = c("solid", "dashed"), 
                        labels = c("Superior   ", "Ontario   ")) +
  scale_y_continuous(limits = c(-2, 60), breaks = seq(0, 70, 10), expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0.1)) +
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
                                     get_legend(plot.m2.es),
                                     nrow = 1,
                                     widths = c(0.09, 1)),
                         arrangeGrob(plot.m2.es + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank()),
                                     plot.m2.dpf + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank()),
                                     plot.m2.add + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank()),
                                     nrow = 3,
                                     left = textGrob("Maternal Effect (%)", y = 0.52, rot = 90, gp = gpar(cex = 1.75, fontfamily = "Arial")),
                                     bottom = textGrob("Incubation Temperature Treatment (°C)", x = 0.545, gp = gpar(cex = 1.75, fontfamily = "Arial"))),
                         heights = c(0.025, 1)
)

ggsave("figures/2020-Embryo-Maternal-SE-Line.png", plot = plot.m2.all, width = 11, height = 15, dpi = 300)

