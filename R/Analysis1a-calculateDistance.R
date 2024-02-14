# Analysis 2a: Assembly processes-calculate distance matrices
# XJ
# 2023.01.29


###########################################################

# clean the workspace

rm(list = ls())

# load library

library(tidyverse)
library(vegan)
library(fields)
library(TeachingDemos)
# source("./R/beta.div.comp.R")
library(adespatial)

char2seed("distances")

# load data
df.soil <- read.csv("./data/data_processing/tarim_soil_data_07212020.csv")
df.xy <- read.csv("./data/data_processing/tarim_xy_locations_12302021.csv")
df.bacteria <- read.csv("./data/data_processing/tarim_bacteria_07212020.csv")
df.fungi <- read.csv("./data/data_processing/tarim_fungi_07212020.csv")


###########################################################
# clean data
df.soil.clean <- df.soil %>% 
  filter(Land_use == "Plantation") %>% 
  droplevels() %>% 
  dplyr::select(id, wells, depth, reps,
                Clay, Silt, Sand, SM, pH, soil_EC,
                soil_Ca, soil_Mg, soil_K, soil_Na, soil_Cl, soil_SO4, soil_HCO3,
                soil_salin,
                SOC, STN, STP, avaP,
                BG, CB, NAG, ALP) %>% 
  inner_join(., df.xy, by = "wells") %>% 
  mutate(wells = factor(wells),
         depth = factor(depth,
                        levels = c("0-10cm", "10-20cm", "20-40cm", 
                                   "40-60cm", "60-80cm", "80-100cm"),
                        labels = c("5", "15", "30", 
                                   "50", "70", "90")),
         reps = factor(reps)) %>% 
  mutate(depth = as.numeric(as.character(depth))) %>% 
  dplyr::select(-id, reps) %>% 
  group_by(wells, depth) %>% 
  summarise_if(is.numeric, mean) %>% 
  rename(well = wells) %>% 
  mutate(well = factor(well,
                       levels = c("s001", "s004", "s007", "s010", "s013",
                                  "s016", "s019", "s023", "s028", "s034"),
                       labels = c("1", "4", "7", "10", "13",
                                  "16", "19", "23", "28", "34"))) %>% 
  mutate(id = paste(well, depth, sep = "_")) %>% 
  mutate(id = factor(id))

df.bacteria.clean <- df.bacteria %>% 
  filter(!treatment == "desert") %>% 
  droplevels() %>% 
  mutate(depth = factor(depth,
                        levels = c("0-10cm", "10-20cm", "20-40cm", "40-60cm",
                                   "60-80cm", "80-100cm"),
                        labels = c("5", "15", "30", 
                                   "50", "70", "90"))) %>% 
  dplyr::select(-X, -treatment, -newtreatment, -rep) %>% 
  mutate(well = factor(well)) %>% 
  pivot_longer(cols = where(is.numeric),
               names_to = "OTUs", values_to = "abundance") %>% 
  group_by(well, depth, OTUs) %>% 
  summarise(abundance = sum(abundance),
            .groups = "drop") %>% 
  pivot_wider(names_from = "OTUs", values_from = "abundance") %>% 
  mutate(id = paste(well, depth, sep = "_")) %>% 
  mutate(id = factor(id)) %>% 
  dplyr::select(-well, -depth)

df.fungi.clean <- df.fungi %>% 
  filter(!treatment == "desert") %>% 
  droplevels() %>% 
  mutate(depth = factor(depth,
                        levels = c("0-10cm", "10-20cm", "20-40cm", "40-60cm",
                                   "60-80cm", "80-100cm"),
                        labels = c("5", "15", "30", 
                                   "50", "70", "90"))) %>% 
  dplyr::select(-X, -treatment, -newtreatment, -rep) %>% 
  mutate(well = factor(well)) %>% 
  pivot_longer(cols = where(is.numeric),
               names_to = "OTUs", values_to = "abundance") %>% 
  group_by(well, depth, OTUs) %>% 
  summarise(abundance = sum(abundance),
            .groups = "drop") %>% 
  pivot_wider(names_from = "OTUs", values_from = "abundance") %>% 
  mutate(depth = as.numeric(as.character(depth)))%>% 
  mutate(id = paste(well, depth, sep = "_")) %>% 
  mutate(id = factor(id)) %>% 
  dplyr::select(-well, -depth)
rm("df.soil", "df.bacteria", "df.fungi", "df.xy")

# double check whether samples are sorted in the same order
# if 0 then OK
sum(df.soil.clean$id != df.bacteria.clean$id)
sum(df.soil.clean$id != df.fungi.clean$id)
sum(df.bacteria.clean$id != df.fungi.clean$id)


# re-clean soil microbial data
# remove sample information
df.bacteria.clean <- df.bacteria.clean %>% 
  dplyr::select(-id)
# remove OTUs of singletons and zeros
df.bacteria.clean <- df.bacteria.clean[, -which(apply(df.bacteria.clean, 2, sum) %in% c(0, 1))]

# remove sample information
df.fungi.clean <- df.fungi.clean %>% 
  dplyr::select(-id)
# remove OTUs of singletons and zeros
df.fungi.clean <- df.fungi.clean[, -which(apply(df.fungi.clean, 2, sum) %in% c(0, 1))]


###########################################################
# data summary for soil microbial data
# # of samples
dim(df.bacteria.clean)[1]
dim(df.fungi.clean)[1]

# # of OTUs
# total
dim(df.bacteria.clean)[2]
dim(df.fungi.clean)[2]
# per samples
sort(apply(df.bacteria.clean > 0, 1, sum))
sort(apply(df.fungi.clean > 0, 1, sum))
# on average
mean(apply(df.bacteria.clean > 0, 1, sum))
mean(apply(df.fungi.clean > 0, 1, sum))

# # of reads
# total
sum(df.bacteria.clean)
sum(df.fungi.clean)
# per samples
sort(apply(df.bacteria.clean, 1, sum))
sort(apply(df.fungi.clean, 1, sum))
# on average
mean(apply(df.bacteria.clean, 1, sum))
mean(apply(df.fungi.clean, 1, sum))


###########################################################
# calculate distances

# geographic distance

geo.dist <- vegdist(df.soil.clean[c("Lat", "Long")],
                    method = "euclidean")

# depth distance

dep.dist <- vegdist(df.soil.clean$depth, method = "euclidean")

# salinity distance

sal.vars <- c("soil_Ca", "soil_Mg", "soil_K", "soil_Na",
              "soil_Cl", "soil_SO4", "soil_HCO3")
sal.dist <- vegdist(apply(df.soil.clean[sal.vars], 2, scale),
                    method = "euclidean")

# nutrient pool distance

nut.pool <- c("SOC", "STN", "avaP")
nut.dist <- vegdist(apply(df.soil.clean[nut.pool], 2, scale),
                    method = "euclidean")

# environmental distance

env.vars <- c("Clay", "Silt", "Sand", "SM", "pH")
env.dist <- vegdist(apply(df.soil.clean[env.vars], 2, scale),
                    method = "euclidean")

# enzymatic distance

enz.vars <- c("BG", "CB", "NAG", "ALP")
enz.dist <- vegdist(apply(df.soil.clean[enz.vars], 2, scale),
                    method = "euclidean")

# microbial beta-diversity

ab.bac.dist <- beta.div.comp(df.bacteria.clean, coef = "S", quant = TRUE)
ab.fun.dist <-  beta.div.comp(df.fungi.clean, coef = "S", quant = TRUE)

# convert to presence/absence scale (0/1)
df.bacteria.pa <- decostand(df.bacteria.clean, method = "pa")
df.fungi.pa <- decostand(df.fungi.clean, method = "pa")
bac.dist <- beta.div.comp(df.bacteria.pa, coef = "S", quant = FALSE)
fun.dist <- beta.div.comp(df.fungi.pa, coef = "S", quant = FALSE)


# combine all the distance matrices
df.dist <- cbind(data.frame(cbind(geo.dist, dep.dist,
                                  sal.dist, env.dist,
                                  nut.dist, enz.dist)),
                 cbind(bac.dist$repl, bac.dist$rich, bac.dist$D),
                 cbind(ab.bac.dist$repl, ab.bac.dist$rich, ab.bac.dist$D),
                 cbind(fun.dist$repl, fun.dist$rich, fun.dist$D),
                 cbind(ab.fun.dist$repl, ab.fun.dist$rich, ab.fun.dist$D))
names(df.dist) <- c("geo", "dep", "sal", "env", "nut", "enz",
                    "bac.repl", "bac.rich", "bac.sor",
                    "bac.repl.wt", "bac.rich.wt", "bac.sor.wt",
                    "fun.repl", "fun.rich", "fun.sor",
                    "fun.repl.wt", "fun.rich.wt", "fun.sor.wt")
m <- bac.dist$repl
m <- as.matrix(m)
m2 <- data.frame(row = rownames(m)[row(m)[lower.tri(m)]],
                 col = colnames(m)[col(m)[lower.tri(m)]])
df.dist <- cbind(m2, df.dist)
write.csv(df.dist, "./data/data_processing/tarim_distance_matrix.csv")


###########################################################
#                    End of Script                        #
###########################################################