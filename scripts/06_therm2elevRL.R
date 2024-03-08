# Aims:
# 1. Calculate elevational range limit corresponding to thermal range limit

# Date created: 13 July 2023
# Date updated: 14 July 2023

# # LIBRARIES # #
library(tidyverse)
library(dplyr)
library(ggplot2)
library(lattice)

rm(list=ls()) 


# # INPUT FILES # #
load('data/census_spp.RData') #adds relative recruitment to census.tms (microhab_analysis_prep.R)


# # OUTPUT FILES # # 
load('data/thermal2elevRL.RData') #thermal to elev range limit per species (therm2elevRL.R)




####################################################################################################

# # EXPLORE DATA # # 

####################################################################################################

load('data/census_spp.RData') #adds relative recruitment to census.tms (microhab_analysis_prep.R)

# Plot all thermal positions
census.spp %>%
  ggplot(aes(x = elev, y = MAT, size = temprange_pos)) + 
  geom_point() 

# Plot thermal positions at upper elev RL
census.spp %>%
  filter(temprange_pos > -1.05 & temprange_pos < -0.95) %>% 
  ggplot(aes(x = elev, y = MAT, size = temprange_pos, color = species)) + 
  geom_point() 

# Plot thermal positions at RL for all species
dat <- census.spp %>% 
  filter(temprange_pos > -1.05 & temprange_pos < -0.95) 
xyplot(MAT ~ elev | species, data = dat)




####################################################################################################

# # DATAFRAME WITH APPROX ELEV RANGE LIMITS # # 

####################################################################################################

# New dataframe at thermal range limit
elev.rl <- census.spp %>% 
  select(full_species, species, MAT, elev, temprange_pos) %>% 
  filter(temprange_pos > -1.05 & temprange_pos < -0.95) %>% 
  distinct() %>% #keep only 1 datapoint per unique row entry for unbiased averaging (i.e. not dependent on how many recruited)
  group_by(species) %>% 
  mutate(mean.elev.rl = mean(elev))

# Check calculations
summary(elev.rl$mean.elev.rl)

elev.rl %>% 
  group_by(species) %>% 
  summarize(min.rl.elev = min(mean.elev.rl), max.rl.elev = max(mean.elev.rl), #these should match
            min.elev = min(elev), mean.elev = mean(elev), max.elev = max(elev), #mean should match above
            min.mat = min(MAT), mean.mat = mean(MAT), max.mat = max(MAT))

# Save DF
save(elev.rl, file = 'data/thermal2elevRL.RData')




