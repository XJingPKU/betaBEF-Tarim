# Path analysis
# XJ
# 2023.02.01

rm(list = ls())

# load library
library(tidyverse)
library(TeachingDemos)
library(lavaan)
library(cowplot)


# set seed
char2seed("pathanalysis")


# load data
df.dist <- read.csv("./data/data_processing/tarim_distance_matrix.csv")

# clean data
df.dist.clean <- df.dist %>% 
  select(-X, -row, -col) %>% 
  mutate(spxd = ifelse(geo != 0 & dep != 0, "SpxD", 
                       ifelse(geo != 0 & dep == 0, "Space", 
                              "Depth"))) %>%  
  mutate(spxd = factor(spxd))

# extract data
overall <- df.dist.clean
space <- df.dist.clean %>% 
  filter(spxd == "Space") %>% 
  droplevels()
depth <- df.dist.clean %>% 
  filter(spxd == "Depth") %>% 
  droplevels()
spxd <- df.dist.clean %>% 
  filter(spxd == "SpxD") %>% 
  droplevels()

# spxd %>% 
#   select(geo, dep, sal, env, nut) %>% 
#   GGally::ggpairs()


###########################################################
# Conduct the path analysis

# replace sor with repl or rich for path models of replacement or richness difference
# space -------------------------------------------------
sem_space <- '
bac.sor ~ geo + sal + env + nut
fun.sor ~ geo + sal + env + nut
enz ~ bac.sor + fun.sor + sal + env + nut
bac.sor ~~ fun.sor
'

fit_space <- sem(sem_space,
                 se = "boot", 
                 bootstrap = 9999,
                 data = space)
summary(fit_space, 
        standardized = TRUE, 
        fit.measures = FALSE, 
        rsquare = TRUE)


# depth -------------------------------------------------
sem_depth <- '
bac.sor ~ dep + sal + env + nut
fun.sor ~ dep + sal + env + nut
enz ~ bac.sor + fun.sor + sal + env + nut
bac.sor ~~ fun.sor
'

fit_depth <- sem(sem_depth,
                 se = "boot", 
                 bootstrap = 9999,
                 data = depth)

summary(fit_depth, 
        standardized = TRUE, 
        fit.measures = FALSE, 
        rsquare = TRUE)


# spxd -------------------------------------------------
sem_spxd <- '
bac.sor ~ geo + dep + sal + env + nut
fun.sor ~ geo + dep + sal + env + nut
enz ~ bac.sor + fun.sor + sal + env + nut
bac.sor ~~ fun.sor
'

fit_spxd <- sem(sem_spxd,
                se = "boot", 
                bootstrap = 9999,
                data = spxd)
summary(fit_spxd, 
        standardized = TRUE, 
        fit.measures = FALSE, 
        rsquare = TRUE)

saveRDS(fit_space, "./outputs/path_analysis_space_sor.RDS")
saveRDS(fit_depth, "./outputs/path_analysis_depth_sor.RDS")
saveRDS(fit_spxd, "./outputs/path_analysis_spxd_sor.RDS")


###########################################################
#                    End of Script                        #
###########################################################