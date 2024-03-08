# Aims:
# 1. Visualize relationships between seedling and microhabitat data
# 2. Choose best random effects structure
# 3. Variance partitioning to check if var is greater between or within sites
# 4. Choose linear vs. quadratic effects

# Date created: 12 Dec 2022
# Date updated: 23 Feb 2023

# # LIBRARIES # #
library(tidyverse)
library(dplyr)
library(ggplot2)
library(lmerTest)
library(AICcmodavg)
library(car)
library(DHARMa)
library(lattice)


rm(list=ls()) 


# # INPUT FILES # #
load('data/census_spp.RData') #adds relative recruitment to census.tms (microhab_analysis_prep.R)
load('data/census_rep.RData') #replicate-level census (2018-2020) and all soil (2022) data (microhab_analysis_prep.R)


# # SUMMARY OF FINDINGS FROM THIS SCRIPT

# MICROHABITAT
# Closed canopies provide a refugia compared to open canopies when looking at average soil moisture 
# and temperature, and plant-height temperature range
# However, variability is high in these relationships
# There is not much correlation between microhabitat variables and elevation
# As expected, many (but not all) of these microhabitat variables are correlated

# SPECIES PATTERNS
# There are species-specific patterns between a few microhabitat variables (e.g. % Fungal) and total 
# recruitment
# Some species are showing a strong response, others a strong varied response, and others little 
# response to env. gradients => need to analyze species by response type!!

# VARIANCE BETWEEN RANDOM EFFECTS
# Variance is greater between sites than between reps
# <1% of the variation attributable to differences between reps within each site (χ = 2.3, p = 0.128) 
# as tested with intraclass correlation test
# 2.7% of variance attributable to differences between sites (χ = 19.2, p < 0.001)
# 4.6% of variance attributable to differences between species (χ = 139.3, p < 0.001)

# MODEL STRUCTURE
# Examining different sensible random effects structures with AICc confirms different species responses,
# as the best structure is: (1|site/rep) + (1|species). This also matches what Katie is using.
# Used AIC to check if quadratic effect provide better fit, and either t3_min or canopycont should 
# be included as a quadratic effect




####################################################################################################

# # PREPPING DATA # # 

####################################################################################################

load('data/census_spp.RData')


# Data frame with species germinating in at least 8 plots -- code from KG (recruitment_optima_analysis.R)
census8 <- census.spp %>% 
  group_by(species) %>% 
  filter(species == "LUPLAT" | species == "ABILAS" | species == "ANEOCC" |
           species == "ERIPER" | species == "PICENG" |
           species == "RUBURS" | species == "SORSIT" | species == "TELGRA" |
           species == "TOLMEN" | species == "VACDEL" | species == "VACPAR" |
           species == "MAHNER" | species == "ERILAN") #"MAHAQU" deleted here because microhabitat data missing for 1 rep 



# Check for correlation in predictor variables
# Predictors with correlation > 0.7 shouldn't be used together (Tabachnick & Fidell (1996))
cc <- c('canopycont', 'microbes', 'f_b', 'perc_fung', 'perc_bact', 'whc', 'c', 'n', 'h',
        't1_avg', 't1_range', 't1_max', 't1_min', 't2_avg', 't2_range', 't2_max', 't2_min',
        't3_avg', 't3_range', 't3_max', 't3_min', 'tmoist_avg', 'tmoist_range', 'tmoist_max', 
        'tmoist_min', 'c_n') #columns of predictor variables

cor(census.spp[cc], use="pairwise.complete.obs")

# Predictor variables with corr > 0.7
# microbes ~ f_b, perc_fung, perc_bact
# f_b ~ perc_fung, perc_bact
# whc ~ c
# c ~ h, n
# h ~ n
# t1_avg ~ t1_range, t2_avg, t2_range, t3_avg
# t2_avg ~ t2_range, t3_avg, t3_range
# t2_range ~ t3_avg, t3_range
# t1_max ~ t1_avg, t1_range, t2_avg, t2_range, t2_max, t3_avg
# t1_min ~ t1_avg, t2_avg, t2_max, t2_min, t3_avg, t3_max
# t2_max ~ t1_avg, t2_avg, t2_range, t3_avg, t3_range, t3_max
# t2_min ~ t3_avg, t3_min, 
# tmoist_avg ~ tmoist_max, tmoist_min

# Conclusion: Chose one variable from each group of variables representing a similar microhabitat
# metric (e.g., f_b, perc_fung, perc_bac, microbes; T1 variables; C-H-N variables) that are 
# uncorrelated. None of the variables are correlated with elevation.
# canopy_cont + f_b + whc + c_n + tmoist_max + tmoist_min + t1_max + t1_min + t3_max + t3_min




####################################################################################################

# # VISUALIZING DATA # # 

####################################################################################################

# # PLOT THEME # #

mytheme <-   theme_classic() +
  theme(axis.text.x = element_text(colour = 'black', size = 25), #x axis text size
        axis.text.y = element_text(colour = 'black', size = 25), #y axis text size
        axis.title.x = element_text(size = 28), #x axis label size
        axis.title.y = element_text(size = 28), #x axis label size
        plot.title = element_text(size = 30, #title size
                                  hjust = 0.5), #align title to center
        legend.title = element_text(size = 24), legend.text = element_text(size = 22)) 

  


####################################################################################################
# Are lowest elevation points are shadier than most? Are closed canopies a refugia?

# Conclusions: Closed canopies are a refugia, especially when looking at soil moisture and temperature ranges. 
# No pattern of soil moisture with elevation or canopy cover. Soil and above-ground temp
# decrease strongly with elevation in RP, but are constant or even increasing with elev in MB. Neither
# soil nor above-ground temp show a pattern with canopy cover.
# Soil moisture range, but not soil moisture average, decreases with elevation with no pattern with 
# canopy cover. Soil moisture avg increases with canopy cover.

dat <- census.rep

# Soil moisture ~ elevation by canopy cover

(gplot <- ggplot(dat, aes(elev, tmoist_avg, color = region)) +
  geom_point(aes(size = canopycont)) +
  ylab('Soil Moisture') + xlab('Elevation [m]') + labs(title = '') + 
  geom_smooth(method = 'lm') +
  mytheme
)


# Soil moisture range ~ elevation by canopy cover
# Conclusion: soil moisture range decreases with elevation

(gplot <- ggplot(dat, aes(elev, tmoist_range, color = region)) +
  geom_point(aes(size = canopycont)) +
  ylab('Soil Moisture Range') + xlab('Elevation [m]') + labs(title = '') + 
  geom_smooth(method = 'lm') +
  mytheme
)


# Soil moisture ~ canopy cover
# Conclusion: soil moisture slightly decreases with canopy cover

(gplot <- ggplot(dat, aes(canopycont, tmoist_avg, color = region)) +
  geom_point() +
  ylab('Soil Moisture') + xlab('Canopy Cover [%]') + labs(title = '') + 
  geom_smooth(method = 'lm') +
  mytheme
)


# Soil moisture range ~ canopy cover
# Conclusion: soil moisture range decreases with canopy cover in MB, vice versa in RP

(gplot <- ggplot(dat, aes(canopycont, tmoist_range, color = region)) +
  geom_point() +
  ylab('Soil Moisture Range') + xlab('Canopy Cover [%]') + labs(title = '') + 
  geom_smooth(method = 'lm') +
  mytheme
)

ggsave('outputs/microhab_Dec2022/tmoist-range_canopy.pdf', gplot)


# Soil temp ~ elevation by canopy cover
# Conclusion: higher elevation areas generally cooler in RP, slightly so for MB

(gplot <- ggplot(dat, aes(elev, t1_avg, color = region)) +
  geom_point(aes(size = canopycont)) +
  ylab('Soil Temperature [ºC]') + xlab('Elevation [m]') + labs(title = '') + 
  geom_smooth(method = 'lm') +
  mytheme
)


# Soil temp ~ canopy cover
# Conclusion: decreases with canopy cover

(gplot <- ggplot(dat, aes(canopycont, t1_avg, color = region)) +
  geom_point() +
  ylab('Soil Temperature [ºC]') + xlab('Canopy Cover [%]') + labs(title = '') + 
  geom_smooth(method = 'lm') +
  mytheme
)


# Soil temp range ~ canopy cover
# Conclusion: decreases in MB, increases in RP

(gplot <- ggplot(dat, aes(canopycont, t1_range, color = region)) +
  geom_point() +
  ylab('Soil Temperature Range [ºC]') + xlab('Canopy Cover [%]') + labs(title = '') + 
  geom_smooth(method = 'lm') +
  mytheme)


# Above ground temp ~ elevation by canopy cover
# Conclusion: temp decreases with elev in RP, but actually increases constant with elev in MB

(gplot <- ggplot(dat, aes(elev, t3_avg, color = region)) +
  geom_point(aes(size = canopycont)) +
  ylab('Above-ground Temperature [ºC]') + xlab('Elevation [m]') + labs(title = '') + 
  geom_smooth(method = 'lm') +
  mytheme
)


# Above ground temp ~ canopy cover
# Conclusion: no real pattern

(gplot <- ggplot(dat, aes(canopycont, t3_avg, color = region)) +
  geom_point() +
  ylab('Above-ground Temperature [ºC]') + xlab('Canopy Cover [%]') + labs(title = '') + 
  geom_smooth(method = 'lm') +
  mytheme
  )


# Above ground temp range ~ elevation by canopy cover
# Conclusion: temp range decreases with elev in RP, but actually increases constant with elev in MB

(gplot <- ggplot(dat, aes(elev, t3_range, color = region)) +
  geom_point(aes(size = canopycont)) +
  ylab('Above-ground Temperature Range [ºC]') + xlab('Elevation [m]') + labs(title = '') +
  geom_smooth(method = 'lm') +
  mytheme
)


# Above ground temp range ~ canopy cover
# Conclusion: increases with elev

(gplot <- ggplot(dat, aes(canopycont, t3_range, color = region)) +
  geom_point() +
  ylab('Above-ground Temperature Range [ºC]') + xlab('Canopy Cover [%]') + labs(title = '') + 
  geom_smooth(method = 'lm') +
  mytheme
)


# Is above ground temp correlated with MAT?
# Conclusion: Yes, in RP but negatively correlated in MB

(gplot <- ggplot(dat, aes(MAT, t3_avg, color = region)) +
  geom_point(aes(size = canopycont)) +
  ylab('Above-ground Temperature [ºC]') + xlab('Mean Annual Temperature [ºC]') + labs(title = '') + 
  geom_smooth(method = 'lm') +
  mytheme
)



####################################################################################################
# Recruitment ~ microhabitat by species
dat <- census.spp

# Soil moisture
(gplot <- ggplot(dat, aes(tmoist_avg, new_all, color = species)) +
geom_point() +
  ylab('Recruitment') + xlab('Soil Moisture') + labs(title = '') +
  geom_smooth(method = 'lm') +
  mytheme)


# Soil temp range
(gplot <- ggplot(dat, aes(t1_range, new_all, color = species)) +
    geom_point() +
    ylab('Recruitment') + xlab('Soil Temperature Range [ºC]') + labs(title = '') +
    geom_smooth(method = 'lm') +
    mytheme)


# Air temp range
(gplot <- ggplot(dat, aes(t3_range, new_all, color = species)) +
    geom_point() +
    ylab('Recruitment') + xlab('Plant-height Temperature Range [ºC]') + labs(title = '') +
    geom_smooth(method = 'lm') +
    mytheme)


# C:N
(gplot <- ggplot(dat, aes(c_n, new_all, color = species)) +
    geom_point() +
    ylab('Recruitment') + xlab('C:N') + labs(title = '') +
    geom_smooth(method = 'lm') +
    mytheme)


# % Fungus
(gplot <- ggplot(dat, aes(perc_fung, new_all, color = species)) +
    geom_point() +
    ylab('Recruitment') + xlab('Soil % Fungal Content') + labs(title = '') +
    geom_smooth(method = 'lm') +
    mytheme)


dat <- census.spp %>% 
  filter(new_all > 0)

# Fungal:Bacterial Ratio
xyplot(new_all ~ f_b | species, data = dat, ylim = c(0, 40), 
       xlab = 'Fungal:Bacterial Ratio', ylab = 'Number of Recruits (outliers removed)') 


# Soil moisture
xyplot(new_all ~ tmoist_min | species, data = dat, ylim = c(0, 40), 
       xlab = 'Minimum Soil Moisture [volumetric water content]', ylab = 'Number of Recruits (outliers removed)') 


# Plant-height temp
xyplot(new_all ~ t3_max | species, data = dat, ylim = c(0, 40), 
       xlab = 'Maximum Plant-Height Temp [ºC]', ylab = 'Number of Recruits (outliers removed)') 


####################################################################################################
# Recruitment ~ microhabitat by species with species germinating in 8 plots or more

dat <- census8

# Canopy cover
xyplot(new_all ~ canopycont | species, data = dat, ylim = c(0, 40)) 

# Fungal:Bacterial Ratio
xyplot(new_all ~ f_b | species, data = dat, ylim = c(0, 40)) 

# WHC
xyplot(new_all ~ whc | species, data = dat, ylim = c(0, 40)) 

# C:N Ratio
xyplot(new_all ~ c_n | species, data = dat, ylim = c(0, 40)) 

# Soil moisture range
xyplot(new_all ~ tmoist_range | species, data = dat, ylim = c(0, 40)) 

# Soil moisture max
xyplot(new_all ~ tmoist_max | species, data = dat, ylim = c(0, 40)) 

# Soil moisture min
xyplot(new_all ~ tmoist_min | species, data = dat, ylim = c(0, 40)) 

# T1 range
xyplot(new_all ~ t1_range | species, data = dat, ylim = c(0, 40)) 

# T1 max
xyplot(new_all ~ t1_max | species, data = dat, ylim = c(0, 40)) 

# T1 min
xyplot(new_all ~ t1_min | species, data = dat, ylim = c(0, 40)) 

# T3 range
xyplot(new_all ~ t3_range | species, data = dat, ylim = c(0, 40)) 

# T3 max
xyplot(new_all ~ t3_max | species, data = dat, ylim = c(0, 40)) 

# T3 min
xyplot(new_all ~ t3_min | species, data = dat, ylim = c(0, 40)) 


# # Conclusions: 
# # - try quadratic variables: tmoist_range, whc, t1_range, t3_range
# # - some species show strong response, others barely any
# # - potentially split by species who had range of germination responses?
# # - patterns with max look very similar to patterns with range, with range patterns looking a bit tighter
# # - patterns with min look different (and min is uncorrelated with range)



####################################################################################################

# # IDENTIFY BEST MODEL STRUCTURE # #  

####################################################################################################

dat <- census.spp

# # Model types

# Fixed Effects (FE): c_n; Random Effects (RE): site/rep + species
mod1 <- lmer(new_all ~ c_n + (1|site/rep) + (1|species), data = dat, REML = F) #REML = F because need log-likelihood for AICc table

# Fixed Effects (FE): c_n; Random Effects (RE): site + species
mod2 <- lmer(new_all ~ c_n + (1|site) + (1|species), data = dat, REML = F) 

# Fixed Effects (FE): c_n; Random Effects (RE): site/rep + family
mod3 <- lmer(new_all ~ c_n + (1|site/rep) + (1|family), data = dat, REML = F) 

# Fixed Effects (FE): c_n; Random Effects (RE): site+ family
mod4 <- lmer(new_all ~ c_n + (1|site) + (1|family), data = dat, REML = F) 

# Fixed Effects (FE): c_n; Random Effects (RE): site/rep
mod5 <- lmer(new_all ~ c_n + (1|site/rep), data = dat, REML = F) #REML = F because need log-likelihood for AICc table

# Fixed Effects (FE): c_n; Random Effects (RE): site
mod6 <- lmer(new_all ~ c_n + (1|site), data = dat, REML = F) 

# Fixed Effects (FE): c_n; Random Effects (RE): site1/rep + species
mod7 <- lmer(new_all ~ c_n + (1|site1/rep) + (1|species), data = dat, REML = F)

# Fixed Effects (FE): c_n; Random Effects (RE): site1/rep + family
mod8 <- lmer(new_all ~ c_n + (1|site1/rep) + (1|family), data = dat, REML = F)

# Fixed Effects (FE): c_n; Random Effects (RE): site/rep + species => doesn't converge
mod1 <- lmer(new_all ~ c_n + (1|site/rep) + (c_n|species), data = dat, REML = F) #REML = F because need log-likelihood for AICc table


# # Sort models by AICc
aictab(list(mod1 = mod1, mod2 = mod2, mod3 = mod3, mod4 = mod4, mod5 = mod5, mod6 = mod6, 
            mod7 = mod7, mod8 = mod8)) 

# # Conclusion: Best model structure is (1|site/rep) + (1|species) OR (1|site) + (1|species), 
# # and given hierarchical sampling design should use (1|site/rep) + (1|species)




####################################################################################################

# # IDENTIFY BEST MODEL STRUCTURE: MULTIVARIATE # #  

####################################################################################################

dat <- census.spp

# # Model types

# Fixed Effects (FE): c_n; Random Effects (RE): site/rep + species
mod1 <- lmer(new_all ~ canopycont + f_b + whc + c_n + tmoist_range + t1_range + t3_range + (1|site/rep) + (1|species), data = dat, REML = F) #REML = F because need log-likelihood for AICc table

# Fixed Effects (FE): c_n; Random Effects (RE): site + species
mod2 <- lmer(new_all ~ canopycont + f_b + whc + c_n + tmoist_range + t1_range + t3_range + (1|site) + (1|species), data = dat, REML = F) 

# Fixed Effects (FE): c_n; Random Effects (RE): site/rep + family
mod3 <- lmer(new_all ~ canopycont + f_b + whc + c_n + tmoist_range + t1_range + t3_range + (1|site/rep) + (1|family), data = dat, REML = F) 

# Fixed Effects (FE): c_n; Random Effects (RE): site+ family
mod4 <- lmer(new_all ~ canopycont + f_b + whc + c_n + tmoist_range + t1_range + t3_range + (1|site) + (1|family), data = dat, REML = F) 

# Fixed Effects (FE): c_n; Random Effects (RE): site/rep
mod5 <- lmer(new_all ~ canopycont + f_b + whc + c_n + tmoist_range + t1_range + t3_range + (1|site/rep), data = dat, REML = F) #REML = F because need log-likelihood for AICc table

# Fixed Effects (FE): c_n; Random Effects (RE): site
mod6 <- lmer(new_all ~ canopycont + f_b + whc + c_n + tmoist_range + t1_range + t3_range + (1|site), data = dat, REML = F) 

# Fixed Effects (FE): c_n; Random Effects (RE): site1/rep + species
mod7 <- lmer(new_all ~ canopycont + f_b + whc + c_n + tmoist_range + t1_range + t3_range + (1|site1/rep) + (1|species), data = dat, REML = F)

# Fixed Effects (FE): c_n; Random Effects (RE): site1/rep + family
mod8 <- lmer(new_all ~ canopycont + f_b + whc + c_n + tmoist_range + t1_range + t3_range + (1|site1/rep) + (1|family), data = dat, REML = F)

# # Sort models by AICc
aictab(list(mod1 = mod1, mod2 = mod2, mod3 = mod3, mod4 = mod4, mod5 = mod5, mod6 = mod6, 
            mod7 = mod7, mod8 = mod8)) 

# # Conclusion: Best model structure is still (1|site/rep) + (1|species) OR (1|site) + (1|species), 
# # and given hierarchical sampling design should use (1|site/rep) + (1|species)




####################################################################################################

# # IDENTIFY BEST MODEL STRUCTURE: UPDATED MULTIVARIATE WITH MAIN RESPONDING SPECIES (CENSUS8) # #  

####################################################################################################

dat <- census.spp %>% 
  group_by(species) %>% 
  filter(species == "LUPLAT" | species == "ABILAS" | species == "ANEOCC" |
           species == "ERIPER" | species == "PICENG" |
           species == "RUBURS" | species == "SORSIT" | species == "TELGRA" |
           species == "TOLMEN" | species == "VACDEL" | species == "VACPAR" |
           species == "MAHNER" | species == "ERILAN") #"MAHAQU" deleted here because microhabitat data missing for 1 rep 


# # Model types

# Random Effects (RE): site/rep + species
mod1 <- lmer(new_all ~ canopycont + f_b + whc + c_n + tmoist_max + tmoist_min + t1_max + 
               t1_min + t3_max + t3_min + (1|site/rep) + (1|species), data = dat, REML = F) #REML = F because need log-likelihood for AICc table

# Random Effects (RE): site + species
mod2 <- lmer(new_all ~ canopycont + f_b + whc + c_n + tmoist_max + tmoist_min + t1_max + 
               t1_min + t3_max + t3_min + (1|site) + (1|species), data = dat, REML = F) 

# Random Effects (RE): site/rep + family
mod3 <- lmer(new_all ~ canopycont + f_b + whc + c_n + tmoist_max + tmoist_min + t1_max + 
               t1_min + t3_max + t3_min + (1|site/rep) + (1|family), data = dat, REML = F) 

# Random Effects (RE): site + family
mod4 <- lmer(new_all ~ canopycont + f_b + whc + c_n + tmoist_max + tmoist_min + t1_max + 
               t1_min + t3_max + t3_min + (1|site) + (1|family), data = dat, REML = F) 

# Random Effects (RE): site/rep
mod5 <- lmer(new_all ~ canopycont + f_b + whc + c_n + tmoist_max + tmoist_min + t1_max + 
               t1_min + t3_max + t3_min + (1|site/rep), data = dat, REML = F) #REML = F because need log-likelihood for AICc table

# Random Effects (RE): site
mod6 <- lmer(new_all ~ canopycont + f_b + whc + c_n + tmoist_max + tmoist_min + t1_max + 
               t1_min + t3_max + t3_min + (1|site), data = dat, REML = F) 

# Random Effects (RE): site1/rep + species
mod7 <- lmer(new_all ~ canopycont + f_b + whc + c_n + tmoist_max + tmoist_min + t1_max + 
               t1_min + t3_max + t3_min + (1|site1/rep) + (1|species), data = dat, REML = F)

# Random Effects (RE): site1/rep + family
mod8 <- lmer(new_all ~ canopycont + f_b + whc + c_n + tmoist_max + tmoist_min + t1_max + 
               t1_min + t3_max + t3_min + (1|site1/rep) + (1|family), data = dat, REML = F)

# Random Effects (RE): site1/rep + rand-slopes species -- failed to converge
mod9 <- lmer(new_all ~ canopycont + f_b + whc + c_n + tmoist_max + tmoist_min + t1_max + 
               t1_min + t3_max + t3_min + (1|site1/rep) + 
              (canopycont + f_b + whc + c_n + tmoist_max + tmoist_min + t1_max + 
                 t1_min + t3_max + t3_min | species), data = dat, REML = F)

# Random Effects (RE): site/rep + rand-slopes species -- failed to converge
mod10 <- lmer(new_all ~ canopycont + f_b + whc + c_n + tmoist_max + tmoist_min + t1_max + 
               t1_min + t3_max + t3_min + (1|site/rep) + 
                (canopycont + f_b + whc + c_n + tmoist_max + tmoist_min + t1_max + 
                   t1_min + t3_max + t3_min | species), data = dat, REML = F)

# Random Effects (RE): site1 + rand-slopes species -- failed to converge
mod11 <- lmer(new_all ~ canopycont + f_b + whc + c_n + tmoist_max + tmoist_min + t1_max + 
               t1_min + t3_max + t3_min + (1|site1) + 
                (canopycont + f_b + whc + c_n + tmoist_max + tmoist_min + t1_max + 
                   t1_min + t3_max + t3_min | species), data = dat, REML = F)

# Random Effects (RE): site + rand-slopes species -- failed to converge
mod12 <- lmer(new_all ~ canopycont + f_b + whc + c_n + tmoist_max + tmoist_min + t1_max + 
                t1_min + t3_max + t3_min + (1|site) + 
                (canopycont + f_b + whc + c_n + tmoist_max + tmoist_min + t1_max + 
                   t1_min + t3_max + t3_min | species), data = dat, REML = F)

# # Sort models by AICc
aictab(list(mod1 = mod1, mod2 = mod2, mod3 = mod3, mod4 = mod4, mod5 = mod5, mod6 = mod6, 
            mod7 = mod7, mod8 = mod8)) 

# # Conclusion: Best model structure is still (1|site/rep) + (1|species) OR (1|site) + (1|species), 
# # and given hierarchical sampling design should use (1|site/rep) + (1|species). 
# # Species random slopes models don't converge => need to define meaningful groupings.




####################################################################################################

# #  INCORPORATE SPECIES RESPONSE INTO MODEL STRUCTURE # #  

####################################################################################################

# Data
dat <- census.spp


# Best model from above - Random Effects (RE): site/rep + species
mod1 <- lmer(new_all ~ canopycont + f_b + whc + c_n + tmoist_max + tmoist_min + t1_max + 
               t1_min + t3_max + t3_min + (1|site/rep) + (1|species), data = dat, REML = F) 

# Model with species as fixed and random effect => boundary(singular) => random effects are too small
mod2 <- lmer(new_all ~ canopycont + f_b + whc + c_n + tmoist_max + tmoist_min + t1_max + 
               t1_min + t3_max + t3_min + species + (1|site/rep) + (1|species), data = dat, REML = F) 

# Model with species as fixed effect only 
mod3 <- lmer(new_all ~ canopycont + f_b + whc + c_n + tmoist_max + tmoist_min + t1_max + 
               t1_min + t3_max + t3_min + species + (1|site/rep), data = dat, REML = F) 

# Model with 'responding' and 'not responding' as fixed effect based off of mod3, and species as RE
spp <- c('ERIPER', 'SORSIT', 'TELGRA', 'TOLMEN', 'VACPAR') #species with slope estimate > 1

dat <- census.spp %>% 
  mutate(respond = if_else(species %in% spp, 1, 0)) #create 'respond' column

mod4 <- lmer(new_all ~ canopycont + f_b + whc + c_n + tmoist_max + tmoist_min + t1_max + 
               t1_min + t3_max + t3_min + respond + (1|site/rep) + (1|species), data = dat, REML = F)

# 'Responding' only as fixed effect and no species RE
mod5 <- lmer(new_all ~ canopycont + f_b + whc + c_n + tmoist_max + tmoist_min + t1_max + 
               t1_min + t3_max + t3_min + respond + (1|site/rep), data = dat, REML = F)

aictab(list(mod1 = mod1, mod2 = mod2, mod3 = mod3, mod4 = mod4, mod5 = mod5))

# # Conclusion: mod 5 is most parsimonious model with 'respond' as fixed effect and no species RE.
# # 'respond' not correlated with any other fixed effects


# Analyze only recruiting species => boundary(singluar) => RE are too small
dat <- census.spp %>% 
  filter(new_all > 0)

mod6 <- lmer(new_all ~ canopycont + f_b + whc + c_n + tmoist_max + tmoist_min + t1_max + 
               t1_min + t3_max + t3_min + (1|site/rep), data = dat, REML = F)

# Analyze only recruiting species with only site as RE => boundary(singluar) => RE are too small
dat <- census.spp %>% 
  filter(new_all > 0)

mod7 <- lmer(new_all ~ canopycont + f_b + whc + c_n + tmoist_max + tmoist_min + t1_max + 
               t1_min + t3_max + t3_min + (1|site), data = dat, REML = F)

# Analyze only recruiting species with only site1 as RE 
dat <- census.spp %>% 
  filter(new_all > 0)

mod8 <- lmer(new_all ~ canopycont + f_b + whc + c_n + tmoist_max + tmoist_min + t1_max + 
               t1_min + t3_max + t3_min + (1|site1), data = dat, REML = F)

# Analyze recruiting vs. not recruiting with only site1 as RE with GLM (family = binomial)
# => model failed to converge
dat <- census.spp %>% 
  mutate(recruit = if_else(new_all > 0, 1, 0)) #create 'recruit' column

mod9 <- glmer(recruit ~ canopycont + f_b + whc + c_n + tmoist_max + tmoist_min + t1_max + 
               t1_min + t3_max + t3_min + (1|site1), data = dat, family = binomial)

#potential solutions: https://stats.stackexchange.com/questions/189115/fitting-a-binomial-glmm-glmer-to-a-response-variable-that-is-a-proportion-or-f




####################################################################################################

# # VARIANCE PARTITIONING # # 
# https://benwhalley.github.io/just-enough-r/icc-and-vpc.html

####################################################################################################

# # Is variance greater between reps or between sites? 
# => Between sites (see below for stats), with majority of variance NOT attributed to RE


# Example model using best model structure identified above

mod <- lmer(new_all ~ c_n + canopycont + t1_avg + f_b + tmoist_avg + 
              (1|site/rep) + (1|species), data = dat) 

# Look at variance estimate

VarCorr(mod)


# Test variance parameters (LRT = chi squared)

rand(mod)


# Calculate variance partition coefficient = variance at a given level of the model, 
# divided by the total variance (the sum of the variance parameters)

VarCorr(mod) %>%
  as_tibble() %>%
  mutate(icc = vcov/sum(vcov)) %>% #variance/total variance
  select(grp, icc)

# => <1% of the variation attributable to differences between reps within each site (χ = 2.3, p = 0.128)
# => 2.7% of variance attributable to differences between sites (χ = 19.2, p < 0.001)
# => 4.6% of variance attributable to differences between species (χ = 139.3, p < 0.001)
# => This means that 91.9% of variance is attributable to Residual variance (not nec. fixed effects!)

# Update 20.3.2023: when tested with site1 (n = 30), variance was higher between replicates than 
# between sites


####################################################################################################

# # CHECK MODEL ASSUMPTIONS (DHARMa code from canopy_analysis.R by KG) # # 
# https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html#general-remarks-on-interperting-residual-patterns-and-tests

####################################################################################################

# Example model based on best model structure identified above

mod <- lmer(new_all ~ c_n + canopycont + t1_avg + f_b + tmoist_avg + 
              (1|site/rep) + (1|species), data = dat) 


# Calculate variance inflation factor (VIF) should be < 4
vif(mod)


# R - Homogeneity of variance: variance in residuals is constant
# L - Assumption of normality: residuals are normally distributed

simres <- simulateResiduals(mod, plot = T)

testUniformity(simres)
testOutliers(simres) 
testQuantiles(simres) 




####################################################################################################

# # IDENTIFY IF SHOULD USE QUADRATIC VS LINEAR EFFECTS # #

####################################################################################################

# Filter dataset to species recruiting in 8 or more plots
dat <- census.spp %>% 
  group_by(species) %>% 
  filter(species == "LUPLAT" | species == "ABILAS" | species == "ANEOCC" |
           species == "ERIPER" | species == "PICENG" |
           species == "RUBURS" | species == "SORSIT" | species == "TELGRA" |
           species == "TOLMEN" | species == "VACDEL" | species == "VACPAR" |
           species == "MAHNER" | species == "ERILAN") #"MAHAQU" deleted here because microhabitat data missing for 1 rep 


# # Fit candidate models

mod1 <- lmer(new_all ~ canopycont + f_b + whc + c_n + tmoist_max + tmoist_min + t1_max + 
               t1_min + t3_max + t3_min + (1|site/rep) + (1|species), data = dat, REML = F) #REML = F because need log-likelihood for AICc table

mod11 <- lmer(new_all ~ canopycont + f_b + whc + c_n + tmoist_max + tmoist_min + t1_max + I(canopycont^2) +
               t1_min + t3_max + t3_min + (1|site/rep) + (1|species), data = dat, REML = F)

mod2 <- lmer(new_all ~ canopycont + f_b + I(f_b^2) + whc + c_n + tmoist_max + tmoist_min + t1_max + 
               t1_min + t3_max + t3_min + (1|site/rep) + (1|species), data = dat, REML = F) 

mod3 <- lmer(new_all ~ canopycont + f_b + whc + I(whc^2) + c_n + tmoist_max + tmoist_min + t1_max + 
               t1_min + t3_max + t3_min + (1|site/rep) + (1|species), data = dat, REML = F)

mod4 <- lmer(new_all ~ canopycont + f_b + whc + c_n + I(c_n^2) + tmoist_max + tmoist_min + t1_max + 
               t1_min + t3_max + t3_min + (1|site/rep) + (1|species), data = dat, REML = F)

mod5 <- lmer(new_all ~ canopycont + f_b + whc + c_n + tmoist_max + I(tmoist_max^2) + tmoist_min + t1_max + 
               t1_min + t3_max + t3_min + (1|site/rep) + (1|species), data = dat, REML = F)

mod6 <- lmer(new_all ~ canopycont + f_b + whc + c_n + tmoist_max + tmoist_min + t1_max + I(tmoist_min^2) +
               t1_min + t3_max + t3_min + (1|site/rep) + (1|species), data = dat, REML = F)

mod7 <- lmer(new_all ~ canopycont + f_b + whc + c_n + tmoist_max + tmoist_min + t1_max + I(t1_min^2) +
               t1_min + t3_max + t3_min + (1|site/rep) + (1|species), data = dat, REML = F)

mod8 <- lmer(new_all ~ canopycont + f_b + whc + c_n + tmoist_max + tmoist_min + t1_max + I(t1_max^2) +
               t1_min + t3_max + t3_min + (1|site/rep) + (1|species), data = dat, REML = F)

mod9 <- lmer(new_all ~ canopycont + f_b + whc + c_n + tmoist_max + tmoist_min + t1_max + I(t3_min^2) +
               t1_min + t3_max + t3_min + (1|site/rep) + (1|species), data = dat, REML = F)

mod10 <- lmer(new_all ~ canopycont + f_b + whc + c_n + tmoist_max + tmoist_min + t1_max + I(t3_max^2) +
               t1_min + t3_max + t3_min + (1|site/rep) + (1|species), data = dat, REML = F)

# # Sort models by AICc
aictab(list(mod1 = mod1, mod2 = mod2, mod3 = mod3, mod4 = mod4, mod5 = mod5, mod6 = mod6, 
            mod7 = mod7, mod8 = mod8, mod9 = mod9, mod10 = mod10, mod11 = mod11)) 


# # Add quadratic effects for canopy and t3_min and compare to baseline model

mod1 <- lmer(new_all ~ canopycont + f_b + whc + c_n + tmoist_max + tmoist_min + t1_max + 
               t1_min + t3_max + t3_min + (1|site/rep) + (1|species), data = dat, REML = F) #REML = F because need log-likelihood for AICc table

mod2 <- lmer(new_all ~ canopycont + f_b + whc + c_n + tmoist_max + tmoist_min + t1_max + 
               I(canopycont^2) + I(t3_min^2) +
                t1_min + t3_max + t3_min + (1|site/rep) + (1|species), data = dat, REML = F)

mod3 <- lmer(new_all ~ canopycont + f_b + whc + c_n + tmoist_max + tmoist_min + t1_max + 
               I(canopycont^2) +
               t1_min + t3_max + t3_min + (1|site/rep) + (1|species), data = dat, REML = F)

mod4 <- lmer(new_all ~ canopycont + f_b + whc + c_n + tmoist_max + tmoist_min + t1_max + 
               I(t3_min^2) +
               t1_min + t3_max + t3_min + (1|site/rep) + (1|species), data = dat, REML = F)

aictab(list(mod1 = mod1, mod2 = mod2, mod3 = mod3, mod4 = mod4))


# # Conclusion: either t3_min or canopycont should be included as a quadratic effect


