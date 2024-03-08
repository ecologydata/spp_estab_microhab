# Aims: MANUSCRIPT TABLES

# Date created: 4 May 2023
# Date updated: 21 Aug 2023


# # LIBRARIES # #
library(tidyverse)
library(dplyr)
library(data.table)


rm(list=ls()) 


# # INPUT FILES # #
load('data/census_spp.RData') #adds relative recruitment to census.tms (microhab_analysis_prep.R)
load('data/census_sdlg.RData') #seedling survival 2018-2022 + microhab data (microhab_analysis-prep.R)
# [rds model files, see below]

load('data/thermal2elevRL.RData') #thermal to elev range limit per species (therm2elevRL.R)

seedmix <- read.csv("data/seedmix_added.csv")


# # OUTPUT FILES # #
load('outputs/mod_avg_bin.RData') #model averaging results for binomial recr (ms_tables.R)
dd.tab <- read.csv('outputs/mod_avg_bin.csv', row.names = F) #model averaging results for binomial recr (ms_tables.R)
load('outputs/mod_avg_bin-cont.RData') #model averaging results for binomial + cont recr (ms_tables.R)
dd.tab <- read.csv('outputs/mod_avg_bin-cont.csv', row.names = F) #model averaging results for binomial + cont recr (ms_tables.R)
load('outputs/ALL_mod_avg.RData') #model averaging results for all models (ms_tables.R)
dd.tab <- read.csv('outputs/ALL_mod_avg.csv') #model averaging results for all models (ms_tables.R)
load('outputs/ALL_mod_avg_SIMPLE.RData') #simplified model avg results to +/- (ms_tables.R)
dd.simple <- read.csv('outputs/ALL_mod_avg_SIMPLE.csv')#simplified model avg results to +/- (ms_tables.R)

elev.limit <- read.csv('outputs/elev-limit.csv') #calculated elev range limit (ms_tables.R)

tab.all.bySpecies <- ('outputs/sites-germ-surv-stats-bySpecies.csv')



####################################################################################################

# # MS STATS # # 

####################################################################################################

# Data
load('data/census_spp.RData') #adds relative recruitment to census.tms (microhab_analysis_prep.R)
load('data/census_sdlg.RData') #seedling survival 2018-2022 + microhab data (microhab_analysis-prep.R)

# Study site characteristics
tab <- census.spp %>% 
  group_by(region) %>% 
  summarise(min.elev = min(elev, na.rm = T), max.elev = max(elev, na.rm = T), #elev
            min.x = min(long, na.rm = T), max.x = max(long, na.rm = T), 
            min.y = min(lat, na.rm = T), max.y = max(lat, na.rm = T), 
            min.mat = min(MAT, na.rm = T), max.mat = max(MAT, na.rm = T), 
            min.map = min(MAP, na.rm = T), max.map = max(MAP, na.rm = T))

# Is seedling survival low because of low number of recruits to start out with (-> demographic stochasticity?)
# -> no pattern of 1- or 2-year seedlings survival with recruitment rate
head(census.sdlg)
plot(surv1yr ~ new_all, data = census.sdlg)
plot(prop18surv20 ~ new_all, data = census.sdlg)





####################################################################################################

# # TABLE 1 [mod-results-simple] # # 

####################################################################################################

# Table with full results
load('outputs/ALL_mod_avg.RData') #model averaging results for all models (ms_tables.R)

# Convert values to +/- signs
dd.simple <- dd.tab

for(i in 1:nrow(dd.simple)) {
  
  for (j in 6:ncol(dd.simple)) { #keep 1-6 as is
    
    if(!is.na(dd.simple[i, j]) & dd.simple[i, j] > 0) { #positive effect sizes
      dd.simple[i, j] <- '+'
    }
    
    else if(!is.na(dd.simple[i, j]) & dd.simple[i, j] < 0) { #negative effect sizes
      dd.simple[i, j] <- '-'
    }
  }
}

# Remove R2
dd.simple <- dd.simple[-c(4:5)]

save(dd.simple, file = 'outputs/ALL_mod_avg_SIMPLE.RData')
write.csv(dd.simple, file = 'outputs/ALL_mod_avg_SIMPLE.csv', row.names = F)




####################################################################################################

# # TABLE S1 [microhabitat] # # 

####################################################################################################

# [Created manually from field guide info and lit reviews in GDocs]




####################################################################################################

# # TABLE S2 [germ-surv] # # 

####################################################################################################

# Data
load('data/census_sdlg.RData') #seedling survival 2018-2022 + microhab data (microhab_analysis-prep.R)

# number of seeds added


# number of germinants (total)
tab <- as.data.frame(setDT(census.sdlg)[, list(sum(new_all, na.rm = T)), by = list(site1, species)])
tab <- tab %>% filter(V1 > 0) #keep only sites where sum(new_all) > 0

tab.sum.bySpecies <- as.data.frame(setDT(tab)[, list(sum(V1)), by = list(species)]) #total number of germinants/species
tab.sum.bySpecies <- tab.sum.bySpecies[order(tab.sum.bySpecies$species),] #order by species alphabetically
colnames(tab.sum.bySpecies) <- c('Species', 'Number of Germinants')

# note: new_all == germinants

# number of survivors yr 1
# 2018-19 survivors
tab <- as.data.frame(setDT(census.sdlg)[, list(sum(surv18to19, na.rm = T)), by = list(site1, species)])
tab <- tab %>% filter(V1 > 0) #keep only sites where sum(surv18to19) > 0

tab.18_19.bySpecies <- as.data.frame(setDT(tab)[, list(sum(V1)), by = list(species)]) #total # y1 survivors/species
#tab.18_19.bySpecies <- tab.18_19.bySpecies[order(- tab.18_19.bySpecies$V1),] #descending order
tab.18_19.bySpecies <- tab.18_19.bySpecies[order(tab.18_19.bySpecies$species),] #order by species alphabetically
colnames(tab.18_19.bySpecies) <- c('Species', '2018-19 survivors')
tab.18_19.bySpecies


# 2019-20 survivors
tab <- as.data.frame(setDT(census.sdlg)[, list(sum(surv19to20, na.rm = T)), by = list(site1, species)])
tab <- tab %>% filter(V1 > 0) #keep only sites where sum(surv19to20) > 0

tab.19_20.bySpecies <- as.data.frame(setDT(tab)[, list(sum(V1)), by = list(species)]) #total # y1 survivors/species
tab.19_20.bySpecies <- tab.19_20.bySpecies[order(tab.19_20.bySpecies$species),] #order by species alphabetically
colnames(tab.19_20.bySpecies) <- c('Species', '2019-20 survivors')
tab.19_20.bySpecies

# number of survivors from 2018-2020
tab <- as.data.frame(setDT(census.sdlg)[, list(sum(surv18to20, na.rm = T)), by = list(site1, species)])
tab <- tab %>% filter(V1 > 0) #keep only sites where sum(surv18to20) > 0

tab.18_20.bySpecies <- as.data.frame(setDT(tab)[, list(sum(V1)), by = list(species)]) #total # y1 survivors/species
tab.18_20.bySpecies <- tab.18_20.bySpecies[order(tab.18_20.bySpecies$species),] #order by species alphabetically
colnames(tab.18_20.bySpecies) <- c('Species', '2018-20 survivors')
tab.18_20.bySpecies

# combine all
tab.all.bySpecies <- full_join(tab.sum.bySpecies, tab.18_19.bySpecies, by = 'Species')
tab.all.bySpecies <- full_join(tab.all.bySpecies, tab.19_20.bySpecies, by = 'Species')
tab.all.bySpecies <- full_join(tab.all.bySpecies, tab.18_20.bySpecies, by = 'Species')

# Add full species names back in
spp <- census.sdlg %>% 
  select(species, full_species) %>% 
  rename(Species = species) %>% 
  filter(!is.na(full_species)) %>% 
  distinct() #dataframe with full species names

tab.all.bySpecies <- full_join(tab.all.bySpecies, spp, by = 'Species')

# number of seeds added by species
seedmix <- read.csv("data/raw/seedmix_added.csv")
str(seedmix)
seedmix$species <- as.factor(seedmix$species)

# subset of data added by weight rather than number of seeds
weight <- seedmix[is.na(seedmix$num_add),]
weight <- seedmix[is.na(seedmix$num_add),]
added <- seedmix[is.na(seedmix$weight_add),]


#species range limits for recruitment
load('data/census_spp.RData') #adds relative recruitment to census.tms (microhab_analysis_prep.R)

# New dataframe at thermal range limit between 1 and -1
elev.rl <- census.spp %>% 
  select(full_species, species, MAT, elev, temprange_pos) %>% 
  filter(temprange_pos > -1 & temprange_pos < 1) %>% 
  distinct() %>% #keep only 1 datapoint per unique row entry for unbiased averaging (i.e. not dependent on how many recruited)
  group_by(species) %>% 
  mutate(mean.elev.rl = mean(elev))

# Check calculations
summary(elev.rl$mean.elev.rl)

mean.elevs <- elev.rl %>% 
  group_by(species) %>% 
  summarize(min.elev = min(elev), max.elev = max(elev))

mean.elevs$range <- paste0(mean.elevs$min.elev, "-", mean.elevs$max.elev)
colnames(mean.elevs)[1] <- 'Species'
tab.all.bySpecies <- full_join(tab.all.bySpecies, mean.elevs[,c(1,4)], by = 'Species')

# New dataframe getting max, min, and mean elevations at which recruitment occured
# repeat for elevation range of seeding survival and add new column to MS tables
elev.species <- census.spp %>% 
  select(full_species, species, elev, site1, rel_rec) %>% 
  filter(rel_rec > 0) %>% # only keep cases where recruitment happened 
  distinct() %>% #keep only 1 datapoint per unique row entry for unbiased averaging (i.e. not dependent on how many recruited)
  group_by(species) 

elev.species.ranges <- elev.species %>% 
  group_by(species) %>% 
  summarize(min.elev = min(elev), max.elev = max(elev))

elev.species.ranges$recruitmentRange <- paste0(elev.species.ranges$min.elev, "-", elev.species.ranges$max.elev)

colnames(elev.species.ranges)[1] <- 'Species'
tab.all.bySpecies <- full_join(tab.all.bySpecies, elev.species.ranges[,c(1,4)], by = 'Species')


### seedling survival elevational range
load('data/census_sdlg.RData')

### elevation range where seedlings for each species survived in yr 1 
elev.sdlg.surv1yr <- census.sdlg %>% 
  select(full_species, species, elev, site1, surv1yr) %>% 
  filter(surv1yr > 0) %>% # only keep cases where survival happened 
  distinct() %>% #keep only 1 datapoint per unique row entry for unbiased averaging (i.e. not dependent on how many recruited)
  group_by(species) 

elev.sdlg.surv1yr.ranges <- elev.sdlg.surv1yr %>% 
  group_by(species) %>% 
  summarize(min.elev = min(elev), max.elev = max(elev))

elev.sdlg.surv1yr.ranges$survivalRange <- paste0(elev.sdlg.surv1yr.ranges$min.elev, "-", elev.sdlg.surv1yr.ranges$max.elev)

colnames(elev.sdlg.surv1yr.ranges)[1] <- 'Species'
tab.all.bySpecies <- full_join(tab.all.bySpecies, elev.sdlg.surv1yr.ranges[,c(1,4)], by = 'Species')

### elevation range where seedlings for each species had a prob of surviving 2 years


elev.sdlg.prop18surv20 <- census.sdlg %>% 
  select(full_species, species, elev, site1, prop18surv20) %>% 
  filter(prop18surv20 > 0) %>% # only keep cases where survival happened 
  distinct() %>% #keep only 1 datapoint per unique row entry for unbiased averaging (i.e. not dependent on how many recruited)
  group_by(species) 

elev.sdlg.prop18surv20.ranges <- elev.sdlg.prop18surv20 %>% 
  group_by(species) %>% 
  summarize(min.elev = min(elev), max.elev = max(elev))

elev.sdlg.prop18surv20.ranges$survivalprop18surv20Range <- paste0(elev.sdlg.prop18surv20.ranges$min.elev, "-", elev.sdlg.prop18surv20.ranges$max.elev)

colnames(elev.sdlg.prop18surv20.ranges)[1] <- 'Species'
tab.all.bySpecies <- full_join(tab.all.bySpecies, elev.sdlg.prop18surv20.ranges[,c(1,4)], by = 'Species')

tab.all.bySpecies <- tab.all.bySpecies %>% 
  select(- Species) %>% 
  rename(Species = full_species) %>% 
  relocate(Species, .before = Germination)

tab.all.bySpecies

# Save
write.csv(tab.all.bySpecies, file = 'outputs/sites-germ-surv-stats-bySpecies.csv', row.names = F)




####################################################################################################

# # TABLE S3 [mod-results] # # 

####################################################################################################

## RECRUITMENT PROBABILITY ## ------------------------------------------------------------------

## MAHNER ##
avgm <- read_rds('outputs/avgm_bin_mahner.rds')
candmods <- read_rds('outputs/candmods_bin_mahner.rds')

# Range of R2 values from candidate models
r2vals <- round(range(candmods$R^2), 2)
R2 <- paste(r2vals[1], r2vals[2], sep = '-')

# Model averaging with delta AICc <= 2 
modavg <- rbind(round(avgm$coefficients[1, ], 2)) #mod avg results
modavg <- as.data.frame(cbind(species = 'MAHNER', type = 'bin', N = nrow(avgm$x), modavg, R2)) #add species
modavg

# Save output as table
dd.temp <- data.frame(species = NA, type = NA, N = NA, regionRP = NA, canopycont = NA, f_b = NA, 
                      whc = NA, c_n = NA, t1_avgmax = NA, tmoist_avgmin = NA, 
                      winteravgmn = NA, springtotalsnow = NA)

dd.tab <- bind_rows(dd.temp, modavg)
dd.tab <- dd.tab[-1, ] #remove initialization row
dd.tab


## ERILAN ##
avgm <- read_rds('outputs/avgm_bin_erilan.rds')
candmods <- read_rds('outputs/candmods_bin_erilan.rds')

# Range of R2 values from candidate models
r2vals <- round(range(candmods$R^2), 2)
R2 <- paste(r2vals[1], r2vals[2], sep = '-')

# Model averaging with delta AICc <= 2 
modavg <- rbind(round(avgm$coefficients[1, ], 2)) #mod avg results
modavg <- as.data.frame(cbind(species = 'ERILAN', type = 'bin', N = nrow(avgm$x), modavg, R2)) #add species
modavg

# Add output to table from last section
dd.tab <- bind_rows(dd.tab, modavg)
dd.tab


## VACPAR ##
avgm <- read_rds('outputs/avgm_bin_vacpar.rds')
candmods <- read_rds('outputs/candmods_bin_vacpar.rds')

# Range of R2 values from candidate models
r2vals <- round(range(candmods$R^2), 2)
R2 <- paste(r2vals[1], r2vals[2], sep = '-')

# Model averaging with delta AICc <= 2 
modavg <- rbind(round(avgm$coefficients[1, ], 2)) #mod avg results
modavg <- as.data.frame(cbind(species = 'VACPAR', type = 'bin', N = nrow(avgm$x), modavg, R2)) #add species
modavg

# Add output to table from last section
dd.tab <- bind_rows(dd.tab, modavg)
dd.tab


## LUPLAT ##
avgm <- read_rds('outputs/avgm_bin_luplat.rds')
candmods <- read_rds('outputs/candmods_bin_luplat.rds')

# Range of R2 values from candidate models
r2vals <- round(range(candmods$R^2), 2)
R2 <- paste(r2vals[1], r2vals[2], sep = '-')

# Model averaging with delta AICc <= 2 
modavg <- rbind(round(avgm$coefficients[1, ], 2)) #mod avg results
modavg <- as.data.frame(cbind(species = 'LUPLAT', type = 'bin', N = nrow(avgm$x), modavg, R2)) #add species
modavg

# Add output to table from last section
dd.tab <- bind_rows(dd.tab, modavg)
dd.tab


## ABILAS ##
avgm <- read_rds('outputs/avgm_bin_abilas.rds')
candmods <- read_rds('outputs/candmods_bin_abilas.rds')

# Range of R2 values from candidate models
r2vals <- round(range(candmods$R^2), 2)
R2 <- paste(r2vals[1], r2vals[2], sep = '-')

# Model averaging with delta AICc <= 2 
modavg <- rbind(round(avgm$coefficients[1, ], 2)) #mod avg results
modavg <- as.data.frame(cbind(species = 'ABILAS', type = 'bin', N = nrow(avgm$x), modavg, R2)) #add species
modavg

# Add output to table from last section
dd.tab <- bind_rows(dd.tab, modavg)
dd.tab


## ANEOCC ##
avgm <- read_rds('outputs/avgm_bin_anneoc.rds')
candmods <- read_rds('outputs/candmods_bin_anneoc.rds')

# Range of R2 values from candidate models
r2vals <- round(range(candmods$R^2), 2)
R2 <- paste(r2vals[1], r2vals[2], sep = '-')

# Model averaging with delta AICc <= 2 
modavg <- rbind(round(avgm$coefficients[1, ], 2)) #mod avg results
modavg <- as.data.frame(cbind(species = 'ANEOCC', type = 'bin', N = nrow(avgm$x), modavg, R2)) #add species
modavg

# Add output to table from last section
dd.tab <- bind_rows(dd.tab, modavg)
dd.tab


## ERIPER ##
avgm <- read_rds('outputs/avgm_bin_eriper.rds')
candmods <- read_rds('outputs/candmods_bin_eriper.rds')

# Range of R2 values from candidate models
r2vals <- round(range(candmods$R^2), 2)
R2 <- paste(r2vals[1], r2vals[2], sep = '-')

# Model averaging with delta AICc <= 2 
modavg <- rbind(round(avgm$coefficients[1, ], 2)) #mod avg results
modavg <- as.data.frame(cbind(species = 'ERIPER', type = 'bin', N = nrow(avgm$x), modavg, R2)) #add species
modavg

# Add output to table from last section
dd.tab <- bind_rows(dd.tab, modavg)
dd.tab


## PICENG ##
avgm <- read_rds('outputs/avgm_bin_piceng.rds')
candmods <- read_rds('outputs/candmods_bin_piceng.rds')

# Range of R2 values from candidate models
r2vals <- round(range(candmods$R^2), 2)
R2 <- paste(r2vals[1], r2vals[2], sep = '-')

# Model averaging with delta AICc <= 2 
modavg <- rbind(round(avgm$coefficients[1, ], 2)) #mod avg results
modavg <- as.data.frame(cbind(species = 'PICENG', type = 'bin', N = nrow(avgm$x), modavg, R2)) #add species
modavg

# Add output to table from last section
dd.tab <- bind_rows(dd.tab, modavg)
dd.tab


## RUBURS ##
avgm <- read_rds('outputs/avgm_bin_ruburs.rds')
candmods <- read_rds('outputs/candmods_bin_ruburs.rds')

# Range of R2 values from candidate models
r2vals <- round(range(candmods$R^2), 2)
R2 <- paste(r2vals[1], r2vals[2], sep = '-')

# Model averaging with delta AICc <= 2 
modavg <- rbind(round(avgm$coefficients[1, ], 2)) #mod avg results
modavg <- as.data.frame(cbind(species = 'RUBURS', type = 'bin', N = nrow(avgm$x), modavg, R2)) #add species
modavg

# Add output to table from last section
dd.tab <- bind_rows(dd.tab, modavg)
dd.tab


## SORSIT ##
avgm <- read_rds('outputs/avgm_bin_sorsit.rds')
candmods <- read_rds('outputs/candmods_bin_sorsit.rds')

# Range of R2 values from candidate models
r2vals <- round(range(candmods$R^2), 2)
R2 <- paste(r2vals[1], r2vals[2], sep = '-')

# Model averaging with delta AICc <= 2 
modavg <- rbind(round(avgm$coefficients[1, ], 2)) #mod avg results
modavg <- as.data.frame(cbind(species = 'SORSIT', type = 'bin', N = nrow(avgm$x), modavg, R2)) #add species
modavg

# Add output to table from last section
dd.tab <- bind_rows(dd.tab, modavg)
dd.tab


## TELGRA ##
avgm <- read_rds('outputs/avgm_bin_telgra.rds')
candmods <- read_rds('outputs/candmods_bin_telgra.rds')

# Range of R2 values from candidate models
r2vals <- round(range(candmods$R^2), 2)
R2 <- paste(r2vals[1], r2vals[2], sep = '-')

# Model averaging with delta AICc <= 2 
modavg <- rbind(round(avgm$coefficients[1, ], 2)) #mod avg results
modavg <- as.data.frame(cbind(species = 'TELGRA', type = 'bin', N = nrow(avgm$x), modavg, R2)) #add species
modavg

# Add output to table from last section
dd.tab <- bind_rows(dd.tab, modavg)
dd.tab


## TOLMEN ##
avgm <- read_rds('outputs/avgm_bin_tolmen.rds')
candmods <- read_rds('outputs/candmods_bin_tolmen.rds')

# Range of R2 values from candidate models
r2vals <- round(range(candmods$R^2), 2)
R2 <- paste(r2vals[1], r2vals[2], sep = '-')

# Model averaging with delta AICc <= 2 
modavg <- rbind(round(avgm$coefficients[1, ], 2)) #mod avg results
modavg <- as.data.frame(cbind(species = 'TOLMEN', type = 'bin', N = nrow(avgm$x), modavg, R2)) #add species
modavg

# Add output to table from last section
dd.tab <- bind_rows(dd.tab, modavg)
dd.tab


## VACDEL ##
avgm <- read_rds('outputs/avgm_bin_vacdel.rds')
candmods <- read_rds('outputs/candmods_bin_vacdel.rds')

# Range of R2 values from candidate models
r2vals <- round(range(candmods$R^2), 2)
R2 <- paste(r2vals[1], r2vals[2], sep = '-')

# Model averaging with delta AICc <= 2 
modavg <- rbind(round(avgm$coefficients[1, ], 2)) #mod avg results
modavg <- as.data.frame(cbind(species = 'VACDEL', type = 'bin', N = nrow(avgm$x), modavg, R2)) #add species
modavg

# Add output to table from last section
dd.tab <- bind_rows(dd.tab, modavg)
dd.tab


## MAHAQU ##
avgm <- read_rds('outputs/avgm_bin_mahaqu.rds')
candmods <- read_rds('outputs/candmods_bin_mahaqu.rds')

# Range of R2 values from candidate models
r2vals <- round(range(candmods$R^2), 2)
R2 <- paste(r2vals[1], r2vals[2], sep = '-')

# Model averaging with delta AICc <= 2 
modavg <- rbind(round(avgm$coefficients[1, ], 2)) #mod avg results
modavg <- as.data.frame(cbind(species = 'MAHAQU', type = 'bin', N = nrow(avgm$x), modavg, R2)) #add species
modavg

# Add output to table from last section
dd.tab <- bind_rows(dd.tab, modavg)
dd.tab


# INTERMITTENT SAVE
save(dd.tab, file = 'outputs/mod_avg_bin.RData')
write.csv(dd.tab, file = 'outputs/mod_avg_bin.csv', row.names = F)


## RELATIVE NUMBER OF RECRUITS ## ------------------------------------------------------------------

## MAHNER ## 
avgm <- read_rds('outputs/avgm_cont_mahner.rds')
candmods <- read_rds('outputs/candmods_cont_mahner.rds')

# Range of R2 values from candidate models
r2vals <- round(range(candmods$R^2), 2)
R2 <- paste(r2vals[1], r2vals[2], sep = '-')

# Model averaging with delta AICc <= 2 
modavg <- rbind(round(avgm$coefficients[1, ], 2)) #mod avg results
modavg <- as.data.frame(cbind(species = 'MAHNER', type = 'cont', N = nrow(avgm$x), modavg, R2)) #add species
modavg

# Add output to table from last section
dd.tab <- bind_rows(dd.tab, modavg)
dd.tab


## ERILAN ##
avgm <- read_rds('outputs/mod_cont_erilan.rds')
candmods <- read_rds('outputs/candmods_cont_erilan.rds')

# Range of R2 values from candidate models
r2vals <- round(range(candmods$R^2), 2)
R2 <- paste(r2vals[1], r2vals[2], sep = '-')

# Best model
ss <- round(summary(avgm)$coefficients[1], 2)
modavg <- as.data.frame(cbind(species = 'ERILAN', type = 'cont*', N = nobs(avgm), `(Intercept)` = ss[1], R2 = R2))
modavg

# Add output to table from last section
dd.tab <- bind_rows(dd.tab, modavg)
dd.tab


## VACPAR ##
avgm <- read_rds('outputs/avgm_cont_vacpar.rds')
candmods <- read_rds('outputs/candmods_cont_vacpar.rds')

# Range of R2 values from candidate models
r2vals <- round(range(candmods$R^2), 2)
R2 <- paste(r2vals[1], r2vals[2], sep = '-')

# Model averaging with delta AICc <= 2 
modavg <- rbind(round(avgm$coefficients[1, ], 2)) #mod avg results
modavg <- as.data.frame(cbind(species = 'VACPAR', type = 'cont', N = nrow(avgm$x), modavg, R2)) #add species
modavg

# Add output to table from last section
dd.tab <- bind_rows(dd.tab, modavg)
dd.tab


## LUPLAT ##
avgm <- read_rds('outputs/avgm_cont_luplat.rds')
candmods <- read_rds('outputs/candmods_cont_luplat.rds')

# Range of R2 values from candidate models
r2vals <- round(range(candmods$R^2), 2)
R2 <- paste(r2vals[1], r2vals[2], sep = '-')

# Model averaging with delta AICc <= 2 
modavg <- rbind(round(avgm$coefficients[1, ], 2)) #mod avg results
modavg <- as.data.frame(cbind(species = 'LUPLAT', type = 'cont', N = nrow(avgm$x), modavg, R2)) #add species
modavg

# Add output to table from last section
dd.tab <- bind_rows(dd.tab, modavg)
dd.tab


## ABILAS ##
avgm <- read_rds('outputs/avgm_cont_abilas.rds')
candmods <- read_rds('outputs/candmods_cont_abilas.rds')

# Range of R2 values from candidate models
r2vals <- round(range(candmods$R^2), 2)
R2 <- paste(r2vals[1], r2vals[2], sep = '-')

# Model averaging with delta AICc <= 2 
modavg <- rbind(round(avgm$coefficients[1, ], 2)) #mod avg results
modavg <- as.data.frame(cbind(species = 'ABILAS', type = 'cont', N = nrow(avgm$x), modavg, R2)) #add species
modavg

# Add output to table from last section
dd.tab <- bind_rows(dd.tab, modavg)
dd.tab


## ANEOCC ##
avgm <- read_rds('outputs/avgm_cont_aneocc.rds')
candmods <- read_rds('outputs/candmods_cont_aneocc.rds')

# Range of R2 values from candidate models
r2vals <- round(range(candmods$R^2), 2)
R2 <- paste(r2vals[1], r2vals[2], sep = '-')

# Model averaging with delta AICc <= 2 
modavg <- rbind(round(avgm$coefficients[1, ], 2)) #mod avg results
modavg <- as.data.frame(cbind(species = 'ANEOCC', type = 'cont', N = nrow(avgm$x), modavg, R2)) #add species
modavg

# Add output to table from last section
dd.tab <- bind_rows(dd.tab, modavg)
dd.tab


## ERIPER ##
avgm <- read_rds('outputs/avgm_cont_eriper.rds')
candmods <- read_rds('outputs/candmods_cont_eriper.rds')

# Range of R2 values from candidate models
r2vals <- round(range(candmods$R^2), 2)
R2 <- paste(r2vals[1], r2vals[2], sep = '-')

# Model averaging with delta AICc <= 2 
modavg <- rbind(round(avgm$coefficients[1, ], 2)) #mod avg results
modavg <- as.data.frame(cbind(species = 'ERIPER', type = 'cont', N = nrow(avgm$x), modavg, R2)) #add species
modavg

# Add output to table from last section
dd.tab <- bind_rows(dd.tab, modavg)
dd.tab


## PICENG ##
avgm <- read_rds('outputs/avgm_cont_piceng.rds')
candmods <- read_rds('outputs/candmods_cont_piceng.rds')

# Range of R2 values from candidate models
r2vals <- round(range(candmods$R^2), 2)
R2 <- paste(r2vals[1], r2vals[2], sep = '-')

# Model averaging with delta AICc <= 2 
modavg <- rbind(round(avgm$coefficients[1, ], 2)) #mod avg results
modavg <- as.data.frame(cbind(species = 'PICENG', type = 'cont', N = nrow(avgm$x), modavg, R2)) #add species
modavg

# Add output to table from last section
dd.tab <- bind_rows(dd.tab, modavg)
dd.tab


## RUBURS ##
avgm <- read_rds('outputs/avgm_cont_ruburs.rds')
candmods <- read_rds('outputs/candmods_cont_ruburs.rds')

# Range of R2 values from candidate models
r2vals <- round(range(candmods$R^2), 2)
R2 <- paste(r2vals[1], r2vals[2], sep = '-')

# Model averaging with delta AICc <= 2 
modavg <- rbind(round(avgm$coefficients[1, ], 2)) #mod avg results
modavg <- as.data.frame(cbind(species = 'RUBURS', type = 'cont', N = nrow(avgm$x), modavg, R2)) #add species
modavg

# Add output to table from last section
dd.tab <- bind_rows(dd.tab, modavg)
dd.tab


## SORSIT ##
avgm <- read_rds('outputs/avgm_cont_sorsit.rds')
candmods <- read_rds('outputs/candmods_cont_sorsit.rds')

# Range of R2 values from candidate models
r2vals <- round(range(candmods$R^2), 2)
R2 <- paste(r2vals[1], r2vals[2], sep = '-')

# Model averaging with delta AICc <= 2 
modavg <- rbind(round(avgm$coefficients[1, ], 2)) #mod avg results
modavg <- as.data.frame(cbind(species = 'SORSIT', type = 'cont', N = nrow(avgm$x), modavg, R2)) #add species
modavg

# Add output to table from last section
dd.tab <- bind_rows(dd.tab, modavg)
dd.tab


## TELGRA ##
avgm <- read_rds('outputs/avgm_cont_telgra.rds')
candmods <- read_rds('outputs/candmods_cont_telgra.rds')

# Range of R2 values from candidate models
r2vals <- round(range(candmods$R^2), 2)
R2 <- paste(r2vals[1], r2vals[2], sep = '-')

# Model averaging with delta AICc <= 2 
modavg <- rbind(round(avgm$coefficients[1, ], 2)) #mod avg results
modavg <- as.data.frame(cbind(species = 'TELGRA', type = 'cont', N = nrow(avgm$x), modavg, R2)) #add species
modavg

# Add output to table from last section
dd.tab <- bind_rows(dd.tab, modavg)
dd.tab


## TOLMEN ##
avgm <- read_rds('outputs/avgm_cont_tolmen.rds')
candmods <- read_rds('outputs/candmods_cont_tolmen.rds')

# Range of R2 values from candidate models
r2vals <- round(range(candmods$R^2), 2)
R2 <- paste(r2vals[1], r2vals[2], sep = '-')

# Model averaging with delta AICc <= 2 
modavg <- rbind(round(avgm$coefficients[1, ], 2)) #mod avg results
modavg <- as.data.frame(cbind(species = 'TOLMEN', type = 'cont', N = nrow(avgm$x), modavg, R2)) #add species
modavg

# Add output to table from last section
dd.tab <- bind_rows(dd.tab, modavg)
dd.tab


## VACDEL ##
avgm <- read_rds('outputs/avgm_cont_vacdel.rds')
candmods <- read_rds('outputs/candmods_cont_vacdel.rds')

# Range of R2 values from candidate models
r2vals <- round(range(candmods$R^2), 2)
R2 <- paste(r2vals[1], r2vals[2], sep = '-')

# Model averaging with delta AICc <= 2 
modavg <- rbind(round(avgm$coefficients[1, ], 2)) #mod avg results
modavg <- as.data.frame(cbind(species = 'VACDEL', type = 'cont', N = nrow(avgm$x), modavg, R2)) #add species
modavg

# Add output to table from last section
dd.tab <- bind_rows(dd.tab, modavg)
dd.tab


## MAHAQU ##
avgm <- read_rds('outputs/avgm_cont_mahaqu.rds')
candmods <- read_rds('outputs/candmods_cont_mahaqu.rds')

# Range of R2 values from candidate models
r2vals <- round(range(candmods$R^2), 2)
R2 <- paste(r2vals[1], r2vals[2], sep = '-')

# Model averaging with delta AICc <= 2 
modavg <- rbind(round(avgm$coefficients[1, ], 2)) #mod avg results
modavg <- as.data.frame(cbind(species = 'MAHAQU', type = 'cont', N = nrow(avgm$x), modavg, R2)) #add species
modavg

# Add output to table from last section
dd.tab <- bind_rows(dd.tab, modavg)
dd.tab


# INTERMITTENT SAVE
save(dd.tab, file = 'outputs/mod_avg_bin-cont.RData')
write.csv(dd.tab, file = 'outputs/mod_avg_bin-cont.csv', row.names = F)


## SEEDLING SURVIVAL ## ------------------------------------------------------------------

## RUBURS ## 
avgm <- read_rds('outputs/avgm_1yr_ruburs.rds')
candmods <- read_rds('outputs/candmods_1yr_ruburs.rds')

# Range of R2 values from candidate models
r2vals <- round(range(candmods$R^2), 2)
R2 <- paste(r2vals[1], r2vals[2], sep = '-')

# Model averaging with delta AICc <= 2 
modavg <- rbind(round(avgm$coefficients[1, ], 2)) #mod avg results
modavg <- as.data.frame(cbind(species = 'RUBURS', ###*** change here
                              type = '1-yr survival', N = nrow(avgm$x), modavg, R2)) 
modavg

# Add output to table from last section
dd.tab <- bind_rows(dd.tab, modavg)
dd.tab[(nrow(dd.tab) - 4):nrow(dd.tab), ]


avgm <- read_rds('outputs/mod_2yr_ruburs.rds')
candmods <- read_rds('outputs/candmods_2yr_ruburs.rds')

# Range of R2 values from candidate models
r2vals <- round(range(candmods$R^2), 2)
R2 <- paste(r2vals[1], r2vals[2], sep = '-')

# Best model summary
ss <- round(summary(avgm)$coefficients[1:4], 2)
modavg <- as.data.frame(cbind(species = 'RUBURS', ###*** change here
                              type = '2-yr survival*', N = nobs(avgm), `(Intercept)` = ss[1],
                        regionRP = ss[2], f_b = ss[3], tmoist_avgmin = ss[4], R2 = R2)) 
modavg

# Add output to table from last section
dd.tab <- bind_rows(dd.tab, modavg)
dd.tab[(nrow(dd.tab) - 4):nrow(dd.tab), ]


## SORSIT ## 
avgm <- read_rds('outputs/avgm_1yr_sorsit.rds')
candmods <- read_rds('outputs/candmods_1yr_sorsit.rds')

# Range of R2 values from candidate models
r2vals <- round(range(candmods$R^2), 2)
R2 <- paste(r2vals[1], r2vals[2], sep = '-')

# Model averaging with delta AICc <= 2 
modavg <- rbind(round(avgm$coefficients[1, ], 2)) #mod avg results
modavg <- as.data.frame(cbind(species = 'SORSIT', ###*** change here
                              type = '1-yr survival', N = nrow(avgm$x), modavg, R2)) 
modavg

# Add output to table from last section
dd.tab <- bind_rows(dd.tab, modavg)
dd.tab[(nrow(dd.tab) - 4):nrow(dd.tab), ]


avgm <- read_rds('outputs/avgm_2yr_sorsit.rds')
candmods <- read_rds('outputs/candmods_2yr_sorsit.rds')

# Range of R2 values from candidate models
r2vals <- round(range(candmods$R^2), 2)
R2 <- paste(r2vals[1], r2vals[2], sep = '-')

# Model averaging with delta AICc <= 2 
modavg <- rbind(round(avgm$coefficients[1, ], 2)) #mod avg results
modavg <- as.data.frame(cbind(species = 'SORSIT', ###*** change here
                              type = '2-yr survival', N = nrow(avgm$x), modavg, R2)) 
modavg

# Add output to table from last section
dd.tab <- bind_rows(dd.tab, modavg)
dd.tab[(nrow(dd.tab) - 4):nrow(dd.tab), ]


## TELGRA ## 
avgm <- read_rds('outputs/avgm_1yr_telgra.rds')
candmods <- read_rds('outputs/candmods_1yr_telgra.rds')

# Range of R2 values from candidate models
r2vals <- round(range(candmods$R^2), 2)
R2 <- paste(r2vals[1], r2vals[2], sep = '-')

# Model averaging with delta AICc <= 2 
modavg <- rbind(round(avgm$coefficients[1, ], 2)) #mod avg results
modavg <- as.data.frame(cbind(species = 'TELGRA', ###*** change here
                              type = '1-yr survival', N = nrow(avgm$x), modavg, R2)) 
modavg

# Add output to table from last section
dd.tab <- bind_rows(dd.tab, modavg)
dd.tab[(nrow(dd.tab) - 4):nrow(dd.tab), ]


avgm <- read_rds('outputs/avgm_2yr_telgra.rds')
candmods <- read_rds('outputs/candmods_2yr_telgra.rds')

# Range of R2 values from candidate models
r2vals <- round(range(candmods$R^2), 2)
R2 <- paste(r2vals[1], r2vals[2], sep = '-')

# Model averaging with delta AICc <= 2 
modavg <- rbind(round(avgm$coefficients[1, ], 2)) #mod avg results
modavg <- as.data.frame(cbind(species = 'TELGRA', ###*** change here
                              type = '2-yr survival', N = nrow(avgm$x), modavg, R2)) 
modavg

# Add output to table from last section
dd.tab <- bind_rows(dd.tab, modavg)
dd.tab[(nrow(dd.tab) - 4):nrow(dd.tab), ]


## TOLMEN ## 
avgm <- read_rds('outputs/avgm_1yr_tolmen.rds')
candmods <- read_rds('outputs/candmods_1yr_tolmen.rds')

# Range of R2 values from candidate models
r2vals <- round(range(candmods$R^2), 2)
R2 <- paste(r2vals[1], r2vals[2], sep = '-')

# Model averaging with delta AICc <= 2 
modavg <- rbind(round(avgm$coefficients[1, ], 2)) #mod avg results
modavg <- as.data.frame(cbind(species = 'TOLMEN', ###*** change here
                              type = '1-yr survival', N = nrow(avgm$x), modavg, R2)) 
modavg

# Add output to table from last section
dd.tab <- bind_rows(dd.tab, modavg)
dd.tab[(nrow(dd.tab) - 4):nrow(dd.tab), ]


avgm <- read_rds('outputs/mod_2yr_tolmen.rds')
candmods <- read_rds('outputs/candmods_2yr_tolmen.rds')

# Range of R2 values from candidate models
r2vals <- round(range(candmods$R^2), 2)
R2 <- paste(r2vals[1], r2vals[2], sep = '-')

# Best model summary
ss <- round(summary(avgm)$coefficients[1:3], 2)
modavg <- as.data.frame(cbind(species = 'TOLMEN', ###*** change here
                              type = '2-yr survival*', N = nobs(avgm), `(Intercept)` = ss[1],
                              c_n = ss[2], f_b = ss[3], R2 = R2)) 
modavg

# Add output to table from last section
dd.tab <- bind_rows(dd.tab, modavg)
dd.tab[(nrow(dd.tab) - 4):nrow(dd.tab), ]


## LUPLAT ## 
avgm <- read_rds('outputs/avgm_1yr_luplat.rds')
candmods <- read_rds('outputs/candmods_1yr_luplat.rds')

# Range of R2 values from candidate models
r2vals <- round(range(candmods$R^2), 2)
R2 <- paste(r2vals[1], r2vals[2], sep = '-')

# Model averaging with delta AICc <= 2 
modavg <- rbind(round(avgm$coefficients[1, ], 2)) #mod avg results
modavg <- as.data.frame(cbind(species = 'LUPLAT', ###*** change here
                              type = '1-yr survival', N = nrow(avgm$x), modavg, R2)) 
modavg

# Add output to table from last section
dd.tab <- bind_rows(dd.tab, modavg)
dd.tab[(nrow(dd.tab) - 4):nrow(dd.tab), ]


avgm <- read_rds('outputs/avgm_2yr_luplat.rds')
candmods <- read_rds('outputs/candmods_2yr_luplat.rds')

# Range of R2 values from candidate models
r2vals <- round(range(candmods$R^2), 2)
R2 <- paste(r2vals[1], r2vals[2], sep = '-')

# Model averaging with delta AICc <= 2 
modavg <- rbind(round(avgm$coefficients[1, ], 2)) #mod avg results
modavg <- as.data.frame(cbind(species = 'LUPLAT', ###*** change here
                              type = '2-yr survival', N = nrow(avgm$x), modavg, R2)) 
modavg

# Add output to table from last section
dd.tab <- bind_rows(dd.tab, modavg)
dd.tab[(nrow(dd.tab) - 4):nrow(dd.tab), ]


## ERIPER ## 
avgm <- read_rds('outputs/avgm_1yr_eriper.rds')
candmods <- read_rds('outputs/candmods_1yr_eriper.rds')

# Range of R2 values from candidate models
r2vals <- round(range(candmods$R^2), 2)
R2 <- paste(r2vals[1], r2vals[2], sep = '-')

# Model averaging with delta AICc <= 2 
modavg <- rbind(round(avgm$coefficients[1, ], 2)) #mod avg results
modavg <- as.data.frame(cbind(species = 'ERIPER', ###*** change here
                              type = '1-yr survival', N = nrow(avgm$x), modavg, R2)) 
modavg

# Add output to table from last section
dd.tab <- bind_rows(dd.tab, modavg)
dd.tab[(nrow(dd.tab) - 4):nrow(dd.tab), ]


avgm <- read_rds('outputs/avgm_2yr_eriper.rds')
candmods <- read_rds('outputs/candmods_2yr_eriper.rds')

# Range of R2 values from candidate models
r2vals <- round(range(candmods$R^2), 2)
R2 <- paste(r2vals[1], r2vals[2], sep = '-')

# Model averaging with delta AICc <= 2 
modavg <- rbind(round(avgm$coefficients[1, ], 2)) #mod avg results
modavg <- as.data.frame(cbind(species = 'ERIPER', ###*** change here
                              type = '2-yr survival', N = nrow(avgm$x), modavg, R2)) 
modavg

# Add output to table from last section
dd.tab <- bind_rows(dd.tab, modavg)
dd.tab[(nrow(dd.tab) - 4):nrow(dd.tab), ]


## VACPAR ## 
avgm <- read_rds('outputs/avgm_1yr_vacpar.rds')
candmods <- read_rds('outputs/candmods_1yr_vacpar.rds')

# Range of R2 values from candidate models
r2vals <- round(range(candmods$R^2), 2)
R2 <- paste(r2vals[1], r2vals[2], sep = '-')

# Model averaging with delta AICc <= 2 
modavg <- rbind(round(avgm$coefficients[1, ], 2)) #mod avg results
modavg <- as.data.frame(cbind(species = 'VACPAR', ###*** change here
                              type = '1-yr survival', N = nrow(avgm$x), modavg, R2)) 
modavg

# Add output to table from last section
dd.tab <- bind_rows(dd.tab, modavg)
dd.tab[(nrow(dd.tab) - 4):nrow(dd.tab), ]


## VACDEL ## 
avgm <- read_rds('outputs/avgm_1yr_vacdel.rds')
candmods <- read_rds('outputs/candmods_1yr_vacdel.rds')

# Range of R2 values from candidate models
r2vals <- round(range(candmods$R^2), 2)
R2 <- paste(r2vals[1], r2vals[2], sep = '-')

# Model averaging with delta AICc <= 2 
modavg <- rbind(round(avgm$coefficients[1, ], 2)) #mod avg results
modavg <- as.data.frame(cbind(species = 'VACDEL', ###*** change here
                              type = '1-yr survival', N = nrow(avgm$x), modavg, R2)) 
modavg

# Add output to table from last section
dd.tab <- bind_rows(dd.tab, modavg)
dd.tab[(nrow(dd.tab) - 4):nrow(dd.tab), ]


## CLEAN UP & SAVE DF ## ------------------------------------------------------------------

load('data/census_spp.RData') #adds relative recruitment to census.tms (microhab_analysis_prep.R)

# Sort by species
dd.tab <- dd.tab[order(dd.tab$species),]

# Add full species names back in
spp <- census.spp %>% 
  select(species, full_species) %>% 
  filter(!is.na(full_species)) %>% 
  distinct() #dataframe with full species names

dd.tab <- left_join(dd.tab, spp, by = 'species')

# Re-arrange columns
dd.tab <- dd.tab[c(28, 2, 3, 16, #model details
                   13, 27, 4, 5, 19, 20, 6, 25, 26, 8, 21, 22, 7, 17, 18, #non-climate params
                   9, 10, 11, 23, 24, 12, 14, 15)]

# Linear estimates of quad effects into one column
dd.tab <- dd.tab %>% 
  mutate(canopycont = if_else(is.na(canopycont), `poly(canopycont, 2)1`, canopycont)) %>% 
  mutate(f_b = if_else(is.na(f_b), `poly(f_b, 2)1`, f_b)) %>% 
  mutate(c_n = if_else(is.na(c_n), `poly(c_n, 2)1`, c_n)) %>% 
  mutate(whc = if_else(is.na(whc), `poly(whc, 2)1`, whc)) %>% 
  mutate(winteravgmn = if_else(is.na(winteravgmn), `poly(winteravgmn, 2)1`, winteravgmn)) %>% 
  mutate(springtotalsnow = if_else(is.na(springtotalsnow), `poly(springtotalsnow, 2)1`, springtotalsnow)) %>% 
  
  #delete duplicate columns
  select(-`poly(canopycont, 2)1`, - `poly(f_b, 2)1`, -`poly(c_n, 2)1`, -`poly(whc, 2)1`, 
         -`poly(winteravgmn, 2)1`, -`poly(springtotalsnow, 2)1`)

# Rename columns
# NOTE: T1 (plant-height) min temp never selected in models so doesn't show up in table
dd.tab <- dd.tab %>%
  rename(Species = full_species, Model = type, Intercept = `(Intercept)`, `Year` = yr1prop19surv20,
         Region = regionRP, `Canopy Openness` = canopycont, `Canopy Openness^2` = `poly(canopycont, 2)2`, 
         `Fungus:Bacteria` = f_b, `Fungus:Bacteria^2` = `poly(f_b, 2)2`, 
         `Carbon:Nitrogen` = c_n, `Carbon:Nitrogen^2` = `poly(c_n, 2)2`,
         `Water Holding Capacity` = whc, `Water Holding Capacity^2` = `poly(whc, 2)2`, 
         `Summer Maximum Soil Temperature` = t1_avgmax, `Summer Minimum Soil Moisture` = tmoist_avgmin,
         `Winter Minimum Soil Temperature` = winteravgmn, `Winter Minimum Soil Temperature^2` = `poly(winteravgmn, 2)2`,
         `Spring Days with Snow` = springtotalsnow, `Spring Days with Snow^2` = `poly(springtotalsnow, 2)2`) %>% 

  mutate(Model = if_else(Model == 'cont', 'number of recruits', Model)) %>% 
  mutate(Model = if_else(Model == 'cont*', 'number of recruits*', Model)) %>% 
  mutate(Model = if_else(Model == 'bin', 'recruitment probability', Model))

# Turn estimates to numeric
dd.tab[5:ncol(dd.tab)] <- apply(dd.tab[5:ncol(dd.tab)], 2, as.numeric)


# FINAL SAVE
save(dd.tab, file = 'outputs/ALL_mod_avg.RData')
write.csv(dd.tab, file = 'outputs/ALL_mod_avg.csv', row.names = F)


