###########################################################
# Supplementary data analysis
# XJ
# 2023.12.14

###########################################################

# clean the workspace 
rm(list = ls())

# load R package 
library(tidyverse)

# load functions
my_se = function(x) {
  se <- sd(x)/sqrt(length(x)) 
}

# load the data 
df.soil <- read.csv("./data/data_processing/tarim_soil_data_07212020.csv")
df.bacteria <- read.csv("./data/data_processing/tarim_bacteria_07212020.csv")
df.bacteria.otu.annotation <- read.csv("./data/data_processing/bacteria_otu_annotation_clean.csv")
df.fungi <- read.csv("./data/data_processing/tarim_fungi_07212020.csv")
df.fungi.otu.annotation <- read.csv("./data/data_processing/fungi_otu_annotation_clean.csv")

###########################################################

# soil characteristics

# data summary
# means
df.soil.clean_mu <- df.soil %>% 
  filter(Land_use == "Plantation") %>% 
  droplevels() %>% 
  dplyr::select(id, wells, depth, reps,
                Clay, Silt, Sand, SM, pH, soil_EC,
                soil_Ca, soil_Mg, soil_K, soil_Na, soil_Cl, soil_SO4, soil_HCO3,
                SOC, STN, avaP,
                BG, CB, NAG, ALP) %>% 
  mutate(soil_sal = soil_Ca + soil_Mg + soil_K + soil_Na +
           soil_Cl + soil_SO4 + soil_HCO3,
         Clay = Clay * 100,
         Silt = Silt * 100,
         Sand = Sand * 100) %>% 
  mutate(wells = factor(wells),
         depth = factor(depth,
                        levels = c("0-10cm", "10-20cm", "20-40cm", 
                                   "40-60cm", "60-80cm", "80-100cm")),
         reps = factor(reps)) %>% 
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
  mutate(id = factor(id)) %>% 
  select(id, well, depth, soil_sal, pH, Clay, Silt, Sand,
         SM, BG, CB, NAG, ALP, SOC, STN, avaP)


# standard errors
df.soil.clean_se <- df.soil %>% 
  filter(Land_use == "Plantation") %>% 
  droplevels() %>% 
  dplyr::select(id, wells, depth, reps,
                Clay, Silt, Sand, SM, pH, soil_EC,
                soil_Ca, soil_Mg, soil_K, soil_Na, soil_Cl, soil_SO4, soil_HCO3,
                SOC, STN, avaP,
                BG, CB, NAG, ALP) %>% 
  mutate(soil_sal = soil_Ca + soil_Mg + soil_K + soil_Na +
           soil_Cl + soil_SO4 + soil_HCO3,
         Clay = Clay * 100,
         Silt = Silt * 100,
         Sand = Sand * 100) %>% 
  mutate(wells = factor(wells),
         depth = factor(depth,
                        levels = c("0-10cm", "10-20cm", "20-40cm", 
                                   "40-60cm", "60-80cm", "80-100cm")),
         reps = factor(reps)) %>% 
  dplyr::select(-id, reps) %>% 
  group_by(wells, depth) %>% 
  summarise_if(is.numeric, my_se) %>% 
  rename(well = wells) %>% 
  mutate(well = factor(well,
                       levels = c("s001", "s004", "s007", "s010", "s013",
                                  "s016", "s019", "s023", "s028", "s034"),
                       labels = c("1", "4", "7", "10", "13",
                                  "16", "19", "23", "28", "34"))) %>% 
  mutate(id = paste(well, depth, sep = "_")) %>% 
  mutate(id = factor(id)) %>% 
  select(id, well, depth, soil_sal, pH, Clay, Silt, Sand,
         SM, BG, CB, NAG, ALP, SOC, STN, avaP)

# save the summary
write.csv(df.soil.clean_mu, "./outputs/soil_data_summary_mu.csv")
write.csv(df.soil.clean_se, "./outputs/soil_data_summary_se.csv")

cor_relation_data <- df.soil.clean_mu %>% 
  select(well, depth, soil_sal, BG, CB, NAG, ALP) %>% 
  pivot_longer(cols = BG:ALP) %>% 
  mutate(name = factor(name,
                       levels = c("BG", "CB", "NAG", "ALP")))

lm_summary <- plyr::ddply(cor_relation_data, 
            c("name", "depth"),
            function(x) {
              temp <- x
              temp$value <- scale(temp$value)
              temp$soil_sal <- scale(temp$soil_sal)
              x <- data.frame(x)
              lmp_mod <- lm(value ~ soil_sal,
                            data = temp)
              suppressWarnings(broom::tidy(lmp_mod)[2, ])
            }
)

write.csv(lm_summary, "./outputs/salinity_enzyme_lm_summary.csv")

lm_summary %>% 
  mutate(depth = forcats::fct_relevel(depth, rev)) %>% 
  ggplot(aes(depth, estimate)) +
  geom_point() +
  geom_linerange(aes(ymin = estimate - 1.96*std.error,
                     ymax = estimate + 1.96*std.error)) +
  geom_hline(yintercept = 0, color = "gray") +
  facet_grid(~ name) +
  coord_flip() +
  labs(x = "Soil depth",
       y = "Standardized regression coefficients") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank())

ggsave("./outputs/salinity_enzyme_lm_summary.pdf",
       width = 180, height = 120, units = "mm")

# cor_relation_data %>% 
#   ggplot(aes(soil_sal, value)) +
#   geom_point() +
#   geom_smooth(method = "lm",
#               formula = "y ~ x",
#               se = FALSE) +
#   facet_grid(name ~ depth, scales = "free_y") +
#   labs(x = expression(paste("Soil salinity (g ", kg^-1, ")")),
#        y = expression(paste("Enzymatic activity (nmol ", g^-1, h^1, ")"))) +
#   theme_bw() +
#   theme(panel.grid = element_blank(),
#         strip.background = element_blank())


###########################################################

# soil microbial community composition
# clean data
df.bacteria.clean <- df.bacteria %>% 
  filter(!treatment == "desert") %>% 
  droplevels() %>% 
  mutate(depth = factor(depth,
                        levels = c("0-10cm", "10-20cm", "20-40cm", 
                                   "40-60cm",
                                   "60-80cm", "80-100cm"))) %>% 
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
  select(id, well, depth, everything())

df.fungi.clean <- df.fungi %>% 
  filter(!treatment == "desert") %>% 
  droplevels() %>% 
  mutate(depth = factor(depth,
                        levels = c("0-10cm", "10-20cm", "20-40cm", "40-60cm",
                                   "60-80cm", "80-100cm"))) %>% 
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
  select(id, well, depth, everything())

# double check whether samples are sorted in the same order
# if 0 then OK
sum(df.soil.clean_mu$id != df.bacteria.clean$id)
sum(df.soil.clean_mu$id != df.fungi.clean$id)
sum(df.bacteria.clean$id != df.fungi.clean$id)

# remove OTUs of singletons and zeros
iid <- df.bacteria.clean %>% 
  select(id, well, depth)
df.bacteria.clean2 <- df.bacteria.clean %>% 
  select(-id, -well, -depth) %>% 
  select(!where(~ sum(.) %in% c(0, 1))) %>% 
  data.frame()

df.fungi.clean2 <- df.fungi.clean %>% 
  select(-id, -well, -depth) %>% 
  select(!where(~ sum(.) %in% c(0, 1))) %>% 
  data.frame()

rownames(df.bacteria.clean2) <- df.soil.clean_mu$id
rownames(df.fungi.clean2) <- df.soil.clean_mu$id

# data summary for soil microbial data
# # of samples
dim(df.bacteria.clean2)[1]
dim(df.fungi.clean2)[1]

# # of OTUs
# total
dim(df.bacteria.clean2)[2]
dim(df.fungi.clean2)[2]
# per samples
sort(apply(df.bacteria.clean2 > 0, 1, sum))
sort(apply(df.fungi.clean2 > 0, 1, sum))
# on average
mean(apply(df.bacteria.clean2 > 0, 1, sum))
mean(apply(df.fungi.clean2 > 0, 1, sum))

# # of reads
# total
sum(df.bacteria.clean2)
sum(df.fungi.clean2)
# per samples
sort(apply(df.bacteria.clean2, 1, sum))
sort(apply(df.fungi.clean2, 1, sum))
# on average
mean(apply(df.bacteria.clean2, 1, sum))
mean(apply(df.fungi.clean2, 1, sum))


df.bacteria.clean2_transpose <- t(df.bacteria.clean2) %>% 
  data.frame() %>% 
  mutate(id = rownames(.))

df.fungi.clean2_transpose <- t(df.fungi.clean2) %>% 
  data.frame()%>% 
  mutate(id = rownames(.))


# linear regressions and composition profiling
# by bacteria
df.bacteria.otu.annotation2 <- df.bacteria.otu.annotation %>% 
  mutate(phylum = paste("a_", phylum)) %>% 
  mutate(phylum = ifelse(phylum == "a_ ", "Unidentified ", phylum)) %>% 
  mutate(phylum = gsub("a_  p__", "", phylum)) %>% 
  mutate(phylum = stringr::str_trim(phylum)) %>% 
  mutate(phylum = factor(phylum)) %>% 
  mutate(phylum = forcats::fct_collapse(phylum,
                                        Acidobacteria = "Acidobacteria",
                                        Actinobacteria = "Actinobacteria",
                                        Bacteroidetes = "Bacteroidetes",
                                        Chloroflexi = "Chloroflexi",
                                        Cyanobacteria = "Cyanobacteria",
                                        Firmicutes = "Firmicutes",
                                        Gemmatimonadetes = "Gemmatimonadetes",
                                        Planctomycetes = "Planctomycetes",
                                        Proteobacteria = "Proteobacteria",
                                        Unidentified  = "Unidentified",
                                        Others = c("[Caldithrix]",
                                                   "[Thermi]",
                                                   "Armatimonadetes",
                                                   "AD3",
                                                   "AncK6",
                                                   "BHI80-139",
                                                   "BRC1",
                                                   "Chlorobi",
                                                   "Chlamydiae",
                                                   "Cyanobacteria",
                                                   "Elusimicrobia",
                                                   "FBP",
                                                   "FCPU426",
                                                   "Fibrobacteres",
                                                   "Fusobacteria",
                                                   "GAL15",
                                                   "GN02",
                                                   "GN04",
                                                   "Kazan-3B-28",
                                                   "Lentisphaerae",
                                                   "MVP-21",
                                                   "NC10",
                                                   "Nitrospirae",
                                                   "NKB19",
                                                   "OD1",
                                                   "OP11",
                                                   "OP3",
                                                   "OP8",
                                                   "PAUC34f",
                                                   "SBR1093",
                                                   "SC4",
                                                   "Spirochaetes",
                                                   "SR1",
                                                   "Tenericutes",
                                                   "TM6",
                                                   "TM7",
                                                   "WPS-2",
                                                   "WS2", 
                                                   "WS3",
                                                   "WWE1",
                                                   "ZB3"))) %>% 
  mutate(phylum = factor(phylum)) %>% 
  select(id, phylum)

  
df.bacteria.clean2_transpose2 <- merge(df.bacteria.clean2_transpose,
                                       df.bacteria.otu.annotation2,
                                       by = "id")

df.bacteria.clean2_phylum <- df.bacteria.clean2_transpose2 %>% 
  group_by(phylum) %>% 
  summarise_if(is.numeric, sum) %>% 
  data.frame()
rownames(df.bacteria.clean2_phylum) <- df.bacteria.clean2_phylum$phylum

df.bacteria.clean2_phylum2 <- df.bacteria.clean2_phylum %>% 
  select(-phylum) %>% 
  t() %>% 
  data.frame()

df.bacteria.clean2_phylum2$site <- rep(c(1, 4, 7, 10, 13, 16, 19, 23, 28, 34), 
                                       each = 6)
df.bacteria.clean2_phylum2$depth <- rep(c("0-10", "10-20", "20-40", "40-60",
                                          "60-80", "80-100"), 10)


df.bacteria.clean2_phylum2_rel_ab <- df.bacteria.clean2_phylum2 %>% 
  mutate(site = factor(site)) %>% 
  group_by(site, depth) %>% 
  summarise_if(is.numeric, sum) %>% 
  pivot_longer(cols = Others:Verrucomicrobia, names_to = "Phylum") %>% 
  mutate(depth = forcats::fct_relevel(depth, rev)) %>% 
  group_by(site, depth) %>% 
  mutate(pct = prop.table(value) * 100) %>% 
  mutate(unique_id = paste(site, depth, sep = "_"))

df.soil.clean_mu2 <- df.soil.clean_mu %>% 
  mutate(depth = gsub("cm", "", depth)) %>% 
  mutate(unique_id = paste(well, depth, sep = "_"))

df.bacteria.clean2_phylum2_rel_ab <- merge(df.bacteria.clean2_phylum2_rel_ab,
                                           df.soil.clean_mu2,
                                           by = "unique_id")

lm_summary_bacteria <- plyr::ddply(df.bacteria.clean2_phylum2_rel_ab, 
                          c("depth.x", "Phylum"),
                          function(x) {
                            temp <- x
                            temp$pct <- scale(temp$pct)
                            temp$soil_sal <- scale(temp$soil_sal)
                            x <- data.frame(x)
                            lmp_mod <- lm(pct ~ soil_sal,
                                          data = temp)
                            suppressWarnings(broom::tidy(lmp_mod)[2, ])
                          }
)
write.csv(lm_summary_bacteria, "./outputs/salinity_bacteria_lm_summary.csv")

lm_summary_bacteria %>% 
  mutate(depth = forcats::fct_relevel(depth.x, rev)) %>% 
  ggplot(aes(depth.x, estimate)) +
  geom_point() +
  geom_linerange(aes(ymin = estimate - 1.96*std.error,
                     ymax = estimate + 1.96*std.error)) +
  geom_hline(yintercept = 0, color = "gray") +
  facet_wrap(~ Phylum) +
  coord_flip() +
  # scale_y_continuous(breaks = c(-1, 0, 1)) +
  labs(x = "Soil depth",
       y = "Standardized regression coefficients") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank())

ggsave("./outputs/salinity_bacteria_lm_summary.pdf",
       width = 180, height = 140, units = "mm")

df.bacteria.clean2_phylum2_rel_ab %>% 
  ggplot(aes(depth.x, pct, color = Phylum, fill = Phylum)) +
  geom_bar(stat = "identity", 
           position = "stack") +
  facet_wrap(~ site, ncol = 3) +
  coord_flip() +
  labs(x = "Soil depth (cm)", 
       y = "Relative Abundance (%)") +
  theme_bw() +
  theme(panel.background = element_blank(),
        strip.background = element_blank())

ggsave("./outputs/salinity_bacteria_profiling.pdf",
       width = 180, height = 140, units = "mm")


# linear regressions and composition profiling
# by fungi
df.fungi.otu.annotation2 <- df.fungi.otu.annotation %>% 
  mutate(phylum = paste("a_", phylum)) %>% 
  mutate(phylum = ifelse(phylum == "a_ ", "Unidentified ", phylum)) %>% 
  mutate(phylum = gsub("a_  p__", "", phylum)) %>% 
  mutate(phylum = stringr::str_trim(phylum)) %>% 
  mutate(phylum = factor(phylum)) %>% 
  mutate(phylum = forcats::fct_collapse(phylum,
                                        Ascomycota = "Ascomycota",
                                        Basidiomycota = "Basidiomycota",
                                        Chytridiomycota = "Chytridiomycota",
                                        Glomeromycota = "Glomeromycota",
                                        Rozellomycota = "Rozellomycota",
                                        Zygomycota = "Zygomycota",
                                        unidentified  = "unidentified",
                                        Others = c("Blastocladiomycota",
                                                   "Neocallimastigomycota"))) %>% 
  mutate(phylum = factor(phylum)) %>% 
  select(id, phylum)


df.fungi.clean2_transpose2 <- merge(df.fungi.clean2_transpose,
                                       df.fungi.otu.annotation2,
                                       by = "id")

df.fungi.clean2_phylum <- df.fungi.clean2_transpose2 %>% 
  group_by(phylum) %>% 
  summarise_if(is.numeric, sum) %>% 
  data.frame()
rownames(df.fungi.clean2_phylum) <- df.fungi.clean2_phylum$phylum

df.fungi.clean2_phylum2 <- df.fungi.clean2_phylum %>% 
  select(-phylum) %>% 
  t() %>% 
  data.frame()

df.fungi.clean2_phylum2$site <- rep(c(1, 4, 7, 10, 13, 16, 19, 23, 28, 34), 
                                       each = 6)
df.fungi.clean2_phylum2$depth <- rep(c("0-10", "10-20", "20-40", "40-60",
                                          "60-80", "80-100"), 10)


df.fungi.clean2_phylum2_rel_ab <- df.fungi.clean2_phylum2 %>% 
  mutate(site = factor(site)) %>% 
  group_by(site, depth) %>% 
  summarise_if(is.numeric, sum) %>% 
  pivot_longer(cols = Ascomycota:Zygomycota, names_to = "Phylum") %>% 
  mutate(depth = forcats::fct_relevel(depth, rev)) %>% 
  group_by(site, depth) %>% 
  mutate(pct = prop.table(value) * 100) %>% 
  mutate(unique_id = paste(site, depth, sep = "_")) %>% 
  filter(Phylum != "Others")


df.fungi.clean2_phylum2_rel_ab <- merge(df.fungi.clean2_phylum2_rel_ab,
                                           df.soil.clean_mu2,
                                           by = "unique_id")

lm_summary_fungi <- plyr::ddply(df.fungi.clean2_phylum2_rel_ab, 
                                   c("depth.x", "Phylum"),
                                   function(x) {
                                     temp <- x
                                     temp$pct <- scale(temp$pct)
                                     temp$soil_sal <- scale(temp$soil_sal)
                                     x <- data.frame(x)
                                     lmp_mod <- lm(pct ~ soil_sal,
                                                   data = temp)
                                     suppressWarnings(broom::tidy(lmp_mod)[2, ])
                                   }
)
write.csv(lm_summary_fungi, "./outputs/salinity_fungi_lm_summary.csv")

lm_summary_fungi %>% 
  # mutate(depth = forcats::fct_relevel(depth.x, rev)) %>% 
  ggplot(aes(depth.x, estimate)) +
  geom_point() +
  geom_linerange(aes(ymin = estimate - 1.96*std.error,
                     ymax = estimate + 1.96*std.error)) +
  geom_hline(yintercept = 0, color = "gray") +
  facet_wrap(~ Phylum) +
  coord_flip() +
  labs(x = "Soil depth",
       y = "Standardized regression coefficients") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank())

ggsave("./outputs/salinity_fungi_lm_summary.pdf",
       width = 180, height = 140, units = "mm")

df.fungi.clean2_phylum2_rel_ab %>% 
  ggplot(aes(depth.x, pct, color = Phylum, fill = Phylum)) +
  geom_bar(stat = "identity", 
           position = "stack") +
  facet_wrap(~ site, ncol = 3) +
  coord_flip() +
  labs(x = "Soil depth (cm)", 
       y = "Relative Abundance (%)") +
  theme_bw() +
  theme(panel.background = element_blank(),
        strip.background = element_blank())

ggsave("./outputs/salinity_fungi_profiling.pdf",
       width = 180, height = 140, units = "mm")


###########################################################
#                    End of Script                        #
###########################################################