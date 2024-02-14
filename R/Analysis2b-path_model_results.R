# Path analysis results
# XJ
# 2023.02.01

rm(list = ls())

# load library
library(tidyverse)
library(lavaan)

# replace sor with repl or rich to get results of replacement or richness differences
# load data
fit_space <- readRDS("./outputs/path_analysis_space_sor.RDS")
fit_depth <- readRDS("./outputs/path_analysis_depth_sor.RDS")
fit_spxd <- readRDS("./outputs/path_analysis_spxd_sor.RDS")

summary(fit_depth, standardized = T, rsq = T, fit.measures = F)
###########################################################
# Summary of the path analysis results

# Standardized path coefficients
nrow_spxd <- c(1:16)
nrow_space <- nrow_depth <- c(1:14)
ncols <- c(1:3, 5:8, 10)

path_coef_std_space <- summary(fit_space, 
                               standardized = TRUE)$PE[nrow_space, ncols] %>% 
  data.frame() %>% 
  mutate(dims = "Horizontal")

path_coef_std_depth <- summary(fit_depth, 
                               standardized = TRUE)$PE[nrow_depth, ncols] %>% 
  data.frame() %>% 
  mutate(dims = "Vertical")

path_coef_std_spxd <- summary(fit_spxd, 
                              standardized = TRUE)$PE[nrow_spxd, ncols] %>% 
  data.frame() %>% 
  mutate(dims = "Horizontal & Vertical")

path_coef_std <- path_coef_std_space %>% 
  bind_rows(., path_coef_std_depth) %>% 
  bind_rows(., path_coef_std_spxd)



# model fit
fit_var <- c("chisq", "df", "pvalue", 
             "cfi", "rmsea", "srmr")
global_fit_measure <- data.frame(rbind(fitMeasures(fit_space, fit_var),
                                       fitMeasures(fit_depth, fit_var),
                                       fitMeasures(fit_spxd, fit_var))) %>% 
  mutate(dims = c("Horizontal",
                  "Vertical",
                  "Horizontal & Vertical"))

# rsquare
mod_rsquare <- data.frame(rbind(cbind(parameterEstimates(fit_space, rsquare = TRUE), 
                                      dims = "Horizontal"),
                                cbind(parameterEstimates(fit_depth, rsquare = TRUE), 
                                      dims = "Vertical"),
                                cbind(parameterEstimates(fit_spxd, rsquare = TRUE), 
                                      dims = "Horizontal & Vertical"))) %>% 
  filter(op == "r2")
mod_rsquare

write.csv(global_fit_measure, "./outputs/global_fit_measures_sor.csv")
write.csv(path_coef_std, "./outputs/path_coef_std_sor.csv")


###########################################################
#                    End of Script                        #
###########################################################