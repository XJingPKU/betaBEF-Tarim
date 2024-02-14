# Analysis 2b: Assembly processes-'varpart' analysis
# XJ
# 2023.01.29


###########################################################

# clean the workspace

rm(list = ls())

# load library

library(tidyverse)
library(TeachingDemos)
library(ecodist)
library(hier.part)
# library(cowplot)

# set seed
char2seed("varpart")

# load functions
source("./R/my_tag_facet.R")

# load data
df.dist <- read.csv("./data/data_processing/tarim_distance_matrix.csv")

# clean data
df.dist.clean <- df.dist %>% 
  select(-X, -row, -col) %>% 
  mutate(spxd = ifelse(geo != 0 & dep != 0, "SpxD", # 60*59/2 = 1770 (1770-270-150 = 1350)
                       ifelse(geo != 0 & dep == 0, "Space", # 10*9/2*6 = 270
                              "Depth"))) %>%  # 6*5/2*10 = 150
  mutate(spxd = factor(spxd))

# extract dimensional data
space <- df.dist.clean %>% 
  filter(spxd == "Space") %>% 
  droplevels()
depth <- df.dist.clean %>% 
  filter(spxd == "Depth") %>% 
  droplevels()
spxd <- df.dist.clean %>% 
  filter(spxd == "SpxD") %>% 
  droplevels()


###########################################################
# inspect salinity distance and beta-diveristy relationships using lmPerm

df_lmp <- df.dist.clean %>% 
  select(spxd, sal, bac.repl, bac.rich, bac.sor,
         fun.repl, fun.rich, fun.sor) %>% 
  pivot_longer(cols = bac.repl:fun.sor, 
               names_to = "variable", 
               values_to = "value") %>% 
  mutate(div_grp = sub("\\.[a-z]+", "", variable),
         div_comp = sub("[a-z]+\\.", "", variable)) %>% 
  mutate(div_grp = factor(div_grp,
                          levels = c("bac", "fun"),
                          labels = c("Bacteria", "Fungi")),
         div_comp = factor(div_comp,
                           levels = c("sor", "repl", "rich")))

res_lmp <- plyr::ddply(df_lmp, 
                       c("spxd", "div_grp", "div_comp"),
                       function(x) {
                         temp <- x
                         temp$value <- scale(temp$value)
                         temp$sal <- scale(temp$sal)
                         x <- data.frame(x)
                         lmp_mod <- lmPerm::lmp(value ~ sal,
                                                data = x,
                                                perm = "Exact",
                                                )
                         suppressWarnings(broom::tidy(lmp_mod)[2, ])
                       }
)


# plot the results
res_reg <- res_lmp %>% 
  mutate(line_type = ifelse(statistic > 0.10, 1, 0),
         estimate = format(round(estimate, 3), nsmall = 3),
         estimate = ifelse(statistic > 0.05 & statistic < 0.10,
                           paste(estimate, ""),
                           # paste(estimate, "\206"), # add dagger symbol by hand
                           ifelse(statistic > 0.10, "--",
                                  estimate))) %>% 
  select(spxd, div_grp, div_comp, estimate, statistic, line_type)

data_text <- res_reg %>% 
  select(spxd, div_grp, div_comp, estimate) %>% 
  mutate(div_scale = factor(div_comp,
                            levels = c("sor", "repl", "rich")),
         x = 8.5,
         y = ifelse(div_comp == "sor", 1.02, 
                    ifelse(div_comp == "repl", .92, .82))) %>% 
  mutate(spxd = factor(spxd,
                       levels = c("Space", "Depth", "SpxD"),
                       labels = c("Horizontal",
                                  "Vertical",
                                  "Horizontal & Vertical")))

reg_plt <- merge(df_lmp, res_reg, by = c("spxd", "div_grp", "div_comp")) %>% 
  mutate(spxd = factor(spxd,
                       levels = c("Space", "Depth", "SpxD"),
                       labels = c("Horizontal",
                                  "Vertical",
                                  "Horizontal & Vertical"))) %>% 
  ggplot(aes(sal, value, color = div_comp)) +
  geom_point(size = .4, alpha = .15) +
  geom_smooth(aes(lty = factor(line_type)),
              method = "lm",
              formula = 'y ~ x',
              se = FALSE,
              lwd = .6) +
  geom_text(aes(x = x, y = y, 
                label = estimate,
                color = div_comp), 
            size = 3.5,
            data = data_text,
            show.legend = FALSE,
            inherit.aes = FALSE) +
  scale_color_manual(values = c("black", "#f46d43", "#74add1"),
                     labels = c(expression(paste(beta["sor"])),
                                expression(paste(beta["repl"])),
                                expression(paste(beta["rich"]))),
                     name = NULL) +
  scale_y_continuous(breaks = c(0, .5, 1), 
                     limits = c(0, 1.05)) +
  scale_linetype_manual(values = c(1, 5),
                        guide = NULL) +
  facet_grid(spxd ~ div_grp) +
  labs(x = "Differences in soil salinity",
       y = expression("Microbial"~beta~"-diversity")) +
  theme_bw(base_size = 11.5) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 11),
        legend.text = element_text(size = 11))
reg_plt
# my_tag_facet(reg_plt)

# save the plot
ggsave("./outputs/betapart_along_salinity_gradients.pdf",
       width = 160, height = 145, units = "mm")


###########################################################
# variance partitioning

# function to extract independent effects of hier.part
extractVarPart <- function(dat, yvar, xvar, dims, type) {
  y <- dat %>% 
    select(all_of(yvar))
  xmat <- dat %>% 
    select(all_of(xvar))
  temp <- hier.part(y, xmat, 
                    gof = "Rsqu", barplot = FALSE)$IJ %>% 
    data.frame()
  I <- temp %>% 
    select(I)
  J <- sum(temp$J)
  res <- rbind(I, covar = J)*100
  res <- cbind(res, dims = dims, type = type, predvar = rownames(res))
  rownames(res) <- NULL
  return(res)
}

# function to extract bootstrapped independent effects of hier.part
boot_extractVarPart <- function(dat, yvar, xvar, dims, type, nboot) {
  varpart_obs <- extractVarPart(dat = dat, 
                                yvar = yvar,
                                xvar = xvar,
                                dims = dims,
                                type = type)
  
  varpart_boot <- dat %>% 
    modelr::bootstrap(n = nboot, id = 'boot_num') %>% 
    group_by(boot_num) %>% 
    summarise(relimp = map(strap, ~ extractVarPart(dat = data.frame(.),
                                                   yvar = yvar,
                                                   xvar = xvar,
                                                   dims = dims,
                                                   type = type))) %>% 
    unnest(relimp) %>% 
    group_by(dims, type, predvar) %>% 
    summarise(lwr_CI = quantile(I, 0.025),
              upr_CI = quantile(I, 0.975),
              .groups = "drop")
  varpart <- merge(varpart_obs, varpart_boot, 
                   by = c("dims", "type", "predvar"))
  return(varpart)
}


# bacteria and fungi --------------------------------------
div_pred_spxd <- c("geo", "dep", "sal", "env", "nut")
div_pred_space <- c("geo", "sal", "env", "nut")
div_pred_depth <- c("dep", "sal", "env", "nut")

# # replace sor with repl or rich for replacement and richness difference
# space_bac_varpart <- boot_extractVarPart(dat = space,
#                                          yvar = "bac.sor",
#                                          xvar = div_pred_space,
#                                          dims = "space",
#                                          type = "bac",
#                                          nboot = 9999)
# depth_bac_varpart <- boot_extractVarPart(dat = depth,
#                                          yvar = "bac.sor",
#                                          xvar = div_pred_depth,
#                                          dims = "depth",
#                                          type = "bac",
#                                          nboot = 9999)
# spxd_bac_varpart <- boot_extractVarPart(dat = spxd,
#                                         yvar = "bac.sor",
#                                         xvar = div_pred_spxd,
#                                         dims = "spxd",
#                                         type = "bac",
#                                         nboot = 9999)
# 
# bac_varpart <- space_bac_varpart %>%
#   bind_rows(., depth_bac_varpart) %>%
#   bind_rows(., spxd_bac_varpart)
# 
# 
# space_fun_varpart <- boot_extractVarPart(dat = space,
#                                          yvar = "fun.sor",
#                                          xvar = div_pred_space,
#                                          dims = "space",
#                                          type = "fun",
#                                          nboot = 9999)
# depth_fun_varpart <- boot_extractVarPart(dat = depth,
#                                          yvar = "fun.sor",
#                                          xvar = div_pred_depth,
#                                          dims = "depth",
#                                          type = "fun",
#                                          nboot = 9999)
# spxd_fun_varpart <- boot_extractVarPart(dat = spxd,
#                                         yvar = "fun.sor",
#                                         xvar = div_pred_spxd,
#                                         dims = "spxd",
#                                         type = "fun",
#                                         nboot = 9999)
# fun_varpart <- space_fun_varpart %>%
#   bind_rows(., depth_fun_varpart) %>%
#   bind_rows(., spxd_fun_varpart)
# 
# 
# div_varpart_clean <- bac_varpart %>%
#   bind_rows(., fun_varpart) %>%
#   mutate(type = factor(type,
#                        levels = c("bac", "fun"),
#                        labels = c("Bacteria", "Fungi")),
#          dims = factor(dims,
#                        levels = c("space", "depth", "spxd"),
#                        labels = c("Horizontal",
#                                   "Vertical", "Horizontal & Vertical")),
#          predvar = factor(predvar,
#                           levels = c("covar", "nut", "env", "sal",
#                                      "dep", "geo")))
# 
# 
# saveRDS(div_varpart_clean, "./outputs/div_varpart_sor.RDS")
# write.csv(div_varpart_clean, "./outputs/BEF_varpart_sor.csv")

# re-load the data
div_varpart_clean <- readRDS("./outputs/div_varpart_sor.RDS")

# plot the results
p_bac <- div_varpart_clean %>% 
  mutate(predvar = factor(predvar,
                          levels = c("covar", "nut", "env", "sal",
                                     "dep", "geo"),
                          labels = c("Covariation",
                                     "Differences in soil nutrient pools",
                                     "Differences in soil properties",
                                     "Differences in soil salinity",
                                     "Differences in soil depth",
                                     "Geographic distance"))) %>% 
  filter(predvar != "Covariation") %>% 
  ggplot(aes(predvar, I)) +
  # geom_col(alpha = .8, fill = "gray70") +
  geom_col(alpha = .3, aes(fill = dims)) +
  geom_pointrange(aes(ymin = lwr_CI, ymax = upr_CI), 
                  shape = 21, fill = "white", stroke = .8) +
  facet_grid(type ~ dims) +
  geom_text(aes(y = 50,
                label = format(round(I, 1), nsmall = 1)),
            size = 3) +
  scale_y_continuous(breaks = seq(0, 50, 20),
                     limits = c(-5, 60)) +
  scale_fill_manual(values = c("#1900b3ff", "#ad1500ff", "#ff9d05ff")) +
  coord_flip() +
  labs(x = NULL, 
       y = "Proportion of explained variance (%)") +
  theme_bw(base_size = 10.5) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        axis.text.y = element_text(size = 9),
        strip.text = element_text(size = 11),
        legend.text = element_text(size = 10),
        legend.position = "none")

my_tag_facet(p_bac, x = 6.5, y = -5)

# save the results
ggsave("./outputs/variance_partitioning_sor.pdf",
       width = 180, height = 140, units = "mm")

###########################################################
#                    End of Script                        #
###########################################################