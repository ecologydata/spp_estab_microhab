# Aims:
# 1. Model averaging with MuMIn for species-level models (1 & 2-yr continuous seedling survival data)

# Date created: 16 June 2023
# Date updated: 12 July 2023

# # LIBRARIES # #
library(tidyverse)
library(dplyr)
library(MuMIn)
library(AICcmodavg)
library(car)
library(data.table)

rm(list=ls()) 


# # INPUT FILES # #
load('data/census_sdlg.RData') #seedling survival 2018-2022 + microhab data (microhab_analysis-prep.R)


# # OUTPUT FILES # #

# Model averaging/best model results
avgm <- read_rds('outputs/avgm_1yr_ruburs.rds') 
candmods <- read_rds('outputs/candmods_1yr_ruburs.rds')

avgm <- read_rds('outputs/mod_2yr_ruburs.rds') #best model
candmods <- read_rds('outputs/candmods_2yr_ruburs.rds')

avgm <- read_rds('outputs/avgm_1yr_sorsit.rds')
candmods <- read_rds('outputs/candmods_1yr_sorsit.rds')

avgm <- read_rds('outputs/avgm_2yr_sorsit.rds')
candmods <- read_rds('outputs/candmods_2yr_sorsit.rds')

avgm <- read_rds('outputs/avgm_1yr_telgra.rds')
candmods <- read_rds('outputs/candmods_1yr_telgra.rds')

avgm <- read_rds('outputs/avgm_2yr_telgra.rds')
candmods <- read_rds('outputs/candmods_2yr_telgra.rds')

avgm <- read_rds('outputs/avgm_1yr_tolmen.rds')
candmods <- read_rds('outputs/candmods_1yr_tolmen.rds')

avgm <- read_rds('outputs/mod_2yr_tolmen.rds') #best model
candmods <- read_rds('outputs/candmods_2yr_tolmen.rds')

avgm <- read_rds('outputs/avgm_1yr_luplat.rds')
candmods <- read_rds('outputs/candmods_1yr_luplat.rds')

avgm <- read_rds('outputs/avgm_2yr_luplat.rds')
candmods <- read_rds('outputs/candmods_2yr_luplat.rds')

avgm <- read_rds('outputs/avgm_1yr_eriper.rds')
candmods <- read_rds('outputs/candmods_1yr_eriper.rds')

avgm <- read_rds('outputs/avgm_2yr_eriper.rds')
candmods <- read_rds('outputs/candmods_2yr_eriper.rds')

avgm <- read_rds('outputs/avgm_1yr_vacpar.rds')
candmods <- read_rds('outputs/candmods_1yr_vacpar.rds')

avgm <- read_rds('outputs/avgm_1yr_vacdel.rds')
candmods <- read_rds('outputs/candmods_1yr_vacdel.rds')




####################################################################################################

# # SET UP DATA FOR MODELS (no TOMST data) # # 

####################################################################################################

# Data
load('data/census_sdlg.RData') #seedling survival 2018-2022 + microhab data (microhab_analysis-prep.R)

# 1-year survival: Species with prop surviving seedlings > 0 per rep in >= 8 plots
tab <- as.data.frame(setDT(census.sdlg)[, list(max(surv1yr, na.rm = T)), by = list(rep, species)])

tab <- tab %>% filter(V1 > 0) #keep only sites where max survival > 0

tab.sum <- as.data.frame(setDT(tab)[, list(length(unique(rep))), by = list(species)]) #unique rep/species
tab.sum <- tab.sum[order(- tab.sum$V1),] #descending order
colnames(tab.sum) <- c('species', 'no.reps')
tab.sum

sdlg.8 <- filter(tab.sum, no.reps >= 8)
sdlg.8.spp <- unique(sdlg.8$species) #SORSIT VACPAR TELGRA TOLMEN RUBURS LUPLAT VACDEL ERIPER

census8.sdlg1yr <- census.sdlg %>% 
  filter(species %in% sdlg.8.spp) 

census8.sdlg1yr$yr1 <- as.factor(census8.sdlg1yr$yr1)


# 2-year survival: Species with prop surviving seedlings > 0 per rep in >= 8 plots
tab <- as.data.frame(setDT(census.sdlg)[, list(max(prop18surv20, na.rm = T)), by = list(rep, species)])

tab <- tab %>% filter(V1 > 0) #keep only sites where max survival > 0

tab.sum <- as.data.frame(setDT(tab)[, list(length(unique(rep))), by = list(species)]) #unique rep/species
tab.sum <- tab.sum[order(- tab.sum$V1),] #descending order
colnames(tab.sum) <- c('species', 'no.reps')
tab.sum

sdlg.8 <- filter(tab.sum, no.reps >= 8)
sdlg.8.spp <- unique(sdlg.8$species) #SORSIT TELGRA TOLMEN ERIPER LUPLAT RUBURS

census8.sdlg2yr <- census.sdlg %>% 
  filter(species %in% sdlg.8.spp)

census8.sdlg2yr$yr2 <- as.factor(census8.sdlg2yr$yr2)


# multi-year survival: max 2 plots where seedling survival > 0
tab <- as.data.frame(setDT(census.sdlg)[, list(max(propxxsurv22, na.rm = T)), by = list(rep, species)])

tab <- tab %>% filter(V1 > 0) #keep only sites where max survival > 0

tab.sum <- as.data.frame(setDT(tab)[, list(length(unique(rep))), by = list(species)]) #unique rep/species
tab.sum <- tab.sum[order(- tab.sum$V1),] #descending order
colnames(tab.sum) <- c('species', 'no.reps')
tab.sum

sdlg.8.spp <- unique(tab.sum$species) 

census8.sdlgxxyr <- census.sdlg %>% 
  filter(species %in% sdlg.8.spp)


# Scale census8.sdlg1yr variables (without TOMST)
census8.sdlg1yr <- census8.sdlg1yr %>%
  
  filter_at(vars(surv1yr, yr1, region, canopycont, f_b, whc, c_n, winteravgmn, springtotalsnow),
            all_vars(!is.na(.))) %>% #remove NAs
  
  mutate_at(vars(canopycont, f_b, whc, c_n, winteravgmn, springtotalsnow), ~scale(.x)[, 1])


# Scale census8.sdlg2yr variables (without TOMST)
census8.sdlg2yr <- census8.sdlg2yr %>%
  
  filter_at(vars(prop18surv20, region, canopycont, f_b, whc, c_n, winteravgmn, springtotalsnow),
            all_vars(!is.na(.))) %>% #remove NAs
  
  mutate_at(vars(canopycont, f_b, whc, c_n, winteravgmn, springtotalsnow), ~scale(.x)[, 1])


# # Full model (with TOMST)
# mod <- lm(surv1yr ~ yr1 + region + canopycont + f_b + whc + c_n + t3_avgmin + tmoist_avgmin + 
#             winteravgmn + springtotalsnow, data = dat) 
# 
# 
# # Full model (without TOMST)
# mod <- lm(surv1yr ~ yr1 + region + canopycont + f_b + whc + c_n + 
#               winteravgmn + springtotalsnow, data = dat) 




####################################################################################################

# # CHECK SEEDLING SURV AT SITES WITH MISSING MICROHAB DATA # # 

####################################################################################################

# Data
load('data/census_sdlg.RData') #seedling survival 2018-2022 + microhab data (microhab_analysis-prep.R)

census8.sdlg1yr <- census.sdlg %>% 
  group_by(species) %>% 
  filter(species == "LUPLAT" | species == "ERIPER" | species == "RUBURS" | species == "SORSIT" | 
           species == "TELGRA" | species == "TOLMEN" | species == "VACDEL" | species == "VACPAR")

# Which variables have many NAs?
census8.sdlg1yr %>% 
  select(canopycont, f_b, whc, c_n, tmoist_avgmin, t1_avgmax, t3_avgmin, winteravgmn, springtotalsnow,
         hobo_avgmax) %>% 
  summary()

# # Conclusion: HOBO > T3 > T1, Tmoist

miss.t3 <- census8.sdlg1yr %>% 
  filter(is.na(t3_avgmin)) %>% #NA values for TOMST T3 sensor
  select(surv1yr, species) %>% 
  group_by(species) %>% 
  summarise(min = min(surv1yr, na.rm = T), mean(surv1yr, na.rm = T), max = max(surv1yr, na.rm = T))

# t3_avgmin NA values don't affect: RUBURS

miss.t1 <- census8.sdlg1yr %>% 
  filter(is.na(t1_avgmax)) %>% #NA values for TOMST T1 sensor (NA for T1 = NA for Tmoist)
  select(surv1yr, species) %>% 
  group_by(species) %>% 
  summarise(min = min(surv1yr, na.rm = T), mean(surv1yr, na.rm = T), max = max(surv1yr, na.rm = T))

# t1_avgmax NA values don't affect: RUBURS

miss.hobo <- census8.sdlg1yr %>% 
  filter(is.na(hobo_avgmax)) %>% #NA values for HOBO
  select(surv1yr, species) %>% 
  group_by(species) %>% 
  summarise(min = min(surv1yr, na.rm = T), mean(surv1yr, na.rm = T), max = max(surv1yr, na.rm = T))

# HOBO NA values don't affect: RUBURS


# # CONCLUSION: run only RUBURS with TOMST variables and remove from all others (3.5.23)




####################################################################################################

# # RUBURS - 1 & 2 YR (with TOMST data) # # *** change here
# Family: Binomial

####################################################################################################

# Data - 1 YEAR SEEDLING SURVIVAL
dat <- census.sdlg %>% 
  
  filter_at(vars(surv1yr, region, canopycont, f_b, whc, c_n, winteravgmn, springtotalsnow, 
                 t3_avgmin, tmoist_avgmin), 
            all_vars(!is.na(.))) %>% #remove NAs
  
  mutate_at(vars(canopycont, f_b, whc, c_n, winteravgmn, springtotalsnow, t3_avgmin, tmoist_avgmin), 
            ~scale(.x)[, 1]) %>% #scale variables
  
  dplyr::filter(species == 'RUBURS')


# Identify how many parameters to use in model selection
mm <- round(nrow(dat)/10, 0) #rule of thumb: ~ 10 data points/parameter
mm

# Choose quadratic effects (<= 2 delta AICc)
mod1 <- glm(surv1yr ~ yr1 + region + poly(canopycont, 2), family = binomial(link = "logit"), data = dat ) #for AICc comparison
mod1.1 <- glm(surv1yr ~ yr1 + region + canopycont, family = binomial(link = "logit") , data = dat )
mod2 <- glm(surv1yr ~ yr1 + region + poly(f_b, 2), family = binomial(link = "logit") , data = dat ) 
mod2.1 <- glm(surv1yr ~ yr1 + region + f_b, family = binomial(link = "logit") , data = dat ) 
mod3 <- glm(surv1yr ~ yr1 + region + poly(whc, 2), family = binomial(link = "logit") , data = dat ) 
mod3.1 <- glm(surv1yr ~ yr1 + region + whc, family = binomial(link = "logit") , data = dat ) 
mod4 <- glm(surv1yr ~ yr1 + region + poly(c_n, 2), family = binomial(link = "logit") , data = dat ) 
mod4.1 <- glm(surv1yr ~ yr1 + region + c_n, family = binomial(link = "logit") , data = dat )
mod5 <- glm(surv1yr ~ yr1 + region + poly(t3_avgmin, 2), family = binomial(link = "logit") , data = dat ) 
mod5.1 <- glm(surv1yr ~ yr1 + region + t3_avgmin, family = binomial(link = "logit") , data = dat ) 
mod7 <- glm(surv1yr ~ yr1 + region + poly(tmoist_avgmin, 2), family = binomial(link = "logit") , data = dat ) 
mod7.1 <- glm(surv1yr ~ yr1 + region + tmoist_avgmin, family = binomial(link = "logit") , data = dat ) 
mod8 <- glm(surv1yr ~ yr1 + region + poly(winteravgmn, 2), family = binomial(link = "logit") , data = dat ) 
mod8.1 <- glm(surv1yr ~ yr1 + region + winteravgmn, family = binomial(link = "logit") , data = dat ) 
mod9 <- glm(surv1yr ~ yr1 + region + poly(springtotalsnow, 2), family = binomial(link = "logit") , data = dat ) 
mod9.1 <- glm(surv1yr ~ yr1 + region + springtotalsnow, family = binomial(link = "logit") , data = dat ) 

aictab(list(canopycont2 = mod1, canopycont = mod1.1, f_b2 = mod2, f_b = mod2.1, whc2 = mod3, whc = mod3.1, 
            c_n2 = mod4, c_n = mod4.1, t1_avgmax2 = mod5, t1_avgmax = mod5.1, 
            tmoist_avgmin2 = mod7, tmoist_avgmin = mod7.1, winteravgmn2 = mod8,
            winteravgmn = mod8.1, springtotalsnow2 = mod9, springtotalsnow = mod9.1))


# Fit global model (= full model)
mod <- glm(surv1yr ~ yr1 + region + canopycont + f_b + whc + c_n + t3_avgmin + tmoist_avgmin + 
             springtotalsnow, family = binomial(link = "logit"), data = dat,
           na.action = "na.fail") 

# Model diagnostics I
vif(mod) 


# # CONCLUSION: removed winteravgmn

# Variable selection with dredge
dd <- dredge(mod, evaluate = T, rank = 'AICc', extra = 'R^2', m.lim = c(0, mm))

# Look at results
subset(dd, delta <= 2)

# Save output to use in microhab_suitab.R
saveRDS(subset(dd, delta <= 2), 'outputs/candmods_1yr_ruburs.rds') 

avgm <- model.avg(get.models(dd, subset = delta <= 2))
saveRDS(avgm, 'outputs/avgm_1yr_ruburs.rds') ###*** change here


# Model diagnostics II
mod <- glm(surv1yr ~ region + canopycont, family = binomial(link = "logit"), data = dat) #fit best model

plot(fitted(mod),resid(mod)) #homogeneity of variance

lag.plot(resid(mod),diag=FALSE,do.lines=FALSE) #independence of residuals to neighboring residuals

deviance(mod)/df.residual(mod) #dispersion should = 1 for variance to equal the mean

coef(lm(with(dat,surv1yr)~fitted(mod))) #link function appropriate if slope = 1

## CONCLUSION: 
# homogeneity of variance - met
# independence of residuals - met
# variance equals the mean - moderately met
# link function - appropriate


# Data - 2 YEAR SEEDLING SURVIVAL
dat <- census.sdlg %>% 
  
  filter_at(vars(prop18surv20, region, canopycont, f_b, whc, c_n, winteravgmn, springtotalsnow, 
                 t3_avgmin, tmoist_avgmin), 
            all_vars(!is.na(.))) %>% #remove NAs
  
  mutate_at(vars(canopycont, f_b, whc, c_n, winteravgmn, springtotalsnow, t3_avgmin, tmoist_avgmin), 
            ~scale(.x)[, 1]) %>% #scale variables
  
  dplyr::filter(species == 'RUBURS')


# Identify how many parameters to use in model selection
mm <- round(nrow(dat)/10, 0) #rule of thumb: ~ 10 data points/parameter
mm

# Choose quadratic effects (<= 2 delta AICc)
mod1 <- glm(prop18surv20 ~ region + poly(canopycont, 2), family = binomial(link = "logit"), data = dat ) #for AICc comparison
mod1.1 <- glm(prop18surv20 ~ region + canopycont, family = binomial(link = "logit") , data = dat )
mod2 <- glm(prop18surv20 ~ region + poly(f_b, 2), family = binomial(link = "logit") , data = dat ) 
mod2.1 <- glm(prop18surv20 ~ region + f_b, family = binomial(link = "logit") , data = dat ) 
mod3 <- glm(prop18surv20 ~ region + poly(whc, 2), family = binomial(link = "logit") , data = dat ) 
mod3.1 <- glm(prop18surv20 ~ region + whc, family = binomial(link = "logit") , data = dat ) 
mod4 <- glm(prop18surv20 ~ region + poly(c_n, 2), family = binomial(link = "logit") , data = dat ) 
mod4.1 <- glm(prop18surv20 ~ region + c_n, family = binomial(link = "logit") , data = dat )
mod5 <- glm(prop18surv20 ~ region + poly(t3_avgmin, 2), family = binomial(link = "logit") , data = dat ) 
mod5.1 <- glm(prop18surv20 ~ region + t3_avgmin, family = binomial(link = "logit") , data = dat ) 
mod7 <- glm(prop18surv20 ~ region + poly(tmoist_avgmin, 2), family = binomial(link = "logit") , data = dat ) 
mod7.1 <- glm(prop18surv20 ~ region + tmoist_avgmin, family = binomial(link = "logit") , data = dat ) 
mod8 <- glm(prop18surv20 ~ region + poly(winteravgmn, 2), family = binomial(link = "logit") , data = dat ) 
mod8.1 <- glm(prop18surv20 ~ region + winteravgmn, family = binomial(link = "logit") , data = dat ) 
mod9 <- glm(prop18surv20 ~ region + poly(springtotalsnow, 2), family = binomial(link = "logit") , data = dat ) 
mod9.1 <- glm(prop18surv20 ~ region + springtotalsnow, family = binomial(link = "logit") , data = dat ) 

aictab(list(canopycont2 = mod1, canopycont = mod1.1, f_b2 = mod2, f_b = mod2.1, whc2 = mod3, whc = mod3.1, 
            c_n2 = mod4, c_n = mod4.1, t1_avgmax2 = mod5, t1_avgmax = mod5.1, 
            tmoist_avgmin2 = mod7, tmoist_avgmin = mod7.1, winteravgmn2 = mod8,
            winteravgmn = mod8.1, springtotalsnow2 = mod9, springtotalsnow = mod9.1))


# Fit global model (= full model)
mod <- glm(prop18surv20 ~ region + f_b + whc + t3_avgmin + tmoist_avgmin, family = binomial(link = "logit"), data = dat,
           na.action = "na.fail") 

# Model diagnostics I
vif(mod) 


# # CONCLUSION: removed poly(tmoist_avgmin), winteravgmn, canopycont, f_b, springtotalsnow, c_n

# Variable selection with dredge
dd <- dredge(mod, evaluate = T, rank = 'AICc', extra = 'R^2', m.lim = c(0, mm))

# Look at results
subset(dd, delta <= 2)

# Model diagnostics II
mod <- glm(prop18surv20 ~ region + f_b + tmoist_avgmin, family = binomial(link = "logit"), 
                  data = dat) #fit best model

# Save best model to use in microhab_suitab.R
saveRDS(subset(dd, delta <= 2), 'outputs/candmods_2yr_ruburs.rds') 

saveRDS(mod, 'outputs/mod_2yr_ruburs.rds') ###*** change here

plot(fitted(mod),resid(mod)) #homogeneity of variance

lag.plot(resid(mod),diag=FALSE,do.lines=FALSE) #independence of residuals to neighboring residuals

deviance(mod)/df.residual(mod) #dispersion should = 1 for variance to equal the mean

coef(lm(with(dat,prop18surv20)~fitted(mod))) #link function appropriate if slope = 1

## CONCLUSION: 
# homogeneity of variance - met
# independence of residuals - met
# variance equals the mean - not met
# link function - appropriate




####################################################################################################

# # SORSIT - 1 & 2 YR # # *** change here
# Family: Binomial

####################################################################################################

# Data - 1 YEAR SEEDLING SURVIVAL
dat <- census8.sdlg1yr %>%
  dplyr::filter(species == 'SORSIT') ###*** change here

# Identify how many parameters to use in model selection
mm <- round(nrow(dat)/10, 0) #rule of thumb: ~ 10 data points/parameter
mm

# Choose quadratic effects (<= 2 delta AICc)
mod1 <- glm(surv1yr ~ yr1 + region + poly(canopycont, 2), family = binomial(link = "logit"), data = dat ) #for AICc comparison
mod1.1 <- glm(surv1yr ~ yr1 + region + canopycont, family = binomial(link = "logit") , data = dat )
mod2 <- glm(surv1yr ~ yr1 + region + poly(f_b, 2), family = binomial(link = "logit") , data = dat ) 
mod2.1 <- glm(surv1yr ~ yr1 + region + f_b, family = binomial(link = "logit") , data = dat ) 
mod3 <- glm(surv1yr ~ yr1 + region + poly(whc, 2), family = binomial(link = "logit") , data = dat ) 
mod3.1 <- glm(surv1yr ~ yr1 + region + whc, family = binomial(link = "logit") , data = dat ) 
mod4 <- glm(surv1yr ~ yr1 + region + poly(c_n, 2), family = binomial(link = "logit") , data = dat ) 
mod4.1 <- glm(surv1yr ~ yr1 + region + c_n, family = binomial(link = "logit") , data = dat )
mod8 <- glm(surv1yr ~ yr1 + region + poly(winteravgmn, 2), family = binomial(link = "logit") , data = dat ) 
mod8.1 <- glm(surv1yr ~ yr1 + region + winteravgmn, family = binomial(link = "logit") , data = dat ) 
mod9 <- glm(surv1yr ~ yr1 + region + poly(springtotalsnow, 2), family = binomial(link = "logit") , data = dat ) 
mod9.1 <- glm(surv1yr ~ yr1 + region + springtotalsnow, family = binomial(link = "logit") , data = dat ) 

aictab(list(canopycont2 = mod1, canopycont = mod1.1, f_b2 = mod2, f_b = mod2.1, whc2 = mod3, whc = mod3.1, 
            c_n2 = mod4, c_n = mod4.1, winteravgmn2 = mod8,
            winteravgmn = mod8.1, springtotalsnow2 = mod9, springtotalsnow = mod9.1))


# Fit global model (= full model)
mod <- glm(surv1yr ~ yr1 + region + poly(canopycont, 2) + f_b + whc + c_n + 
             winteravgmn + poly(springtotalsnow, 2), family = binomial(link = "logit"), data = dat,
           na.action = "na.fail") 

# Model diagnostics I
vif(mod) 

# # CONCLUSION: removed no params

# Variable selection with dredge
dd <- dredge(mod, evaluate = T, rank = 'AICc', extra = 'R^2', m.lim = c(0, mm))

# Look at results
subset(dd, delta <= 2)

# Save output to use in microhab_suitab.R
saveRDS(subset(dd, delta <= 2), 'outputs/candmods_1yr_sorsit.rds') 

avgm <- model.avg(get.models(dd, subset = delta <= 2))
saveRDS(avgm, 'outputs/avgm_1yr_sorsit.rds') ###*** change here


# Model diagnostics II
mod <- glm(surv1yr ~ poly(canopycont, 2) + poly(springtotalsnow, 2) + region,
           family = binomial(link = "logit"), data = dat) #fit best model

plot(fitted(mod),resid(mod)) #homogeneity of variance

lag.plot(resid(mod),diag=FALSE,do.lines=FALSE) #independence of residuals to neighboring residuals

deviance(mod)/df.residual(mod) #dispersion should = 1 for variance to equal the mean

coef(lm(with(dat,surv1yr)~fitted(mod))) #link function appropriate if slope = 1

## CONCLUSION: 
# homogeneity of variance - met
# independence of residuals - met
# variance equals the mean - not met
# link function - appropriate


# Data - 2 YEAR SEEDLING SURVIVAL
dat <- census8.sdlg2yr %>%
  dplyr::filter(species == 'SORSIT') ###*** change here

# Identify how many parameters to use in model selection
mm <- round(nrow(dat)/10, 0) #rule of thumb: ~ 10 data points/parameter
mm

# Choose quadratic effects (<= 2 delta AICc)
mod1 <- glm(prop18surv20 ~ region + poly(canopycont, 2), family = binomial(link = "logit"), data = dat ) #for AICc comparison
mod1.1 <- glm(prop18surv20 ~ region + canopycont, family = binomial(link = "logit") , data = dat )
mod2 <- glm(prop18surv20 ~ region + poly(f_b, 2), family = binomial(link = "logit") , data = dat ) 
mod2.1 <- glm(prop18surv20 ~ region + f_b, family = binomial(link = "logit") , data = dat ) 
mod3 <- glm(prop18surv20 ~ region + poly(whc, 2), family = binomial(link = "logit") , data = dat ) 
mod3.1 <- glm(prop18surv20 ~ region + whc, family = binomial(link = "logit") , data = dat ) 
mod4 <- glm(prop18surv20 ~ region + poly(c_n, 2), family = binomial(link = "logit") , data = dat ) 
mod4.1 <- glm(prop18surv20 ~ region + c_n, family = binomial(link = "logit") , data = dat )
mod8 <- glm(prop18surv20 ~ region + poly(winteravgmn, 2), family = binomial(link = "logit") , data = dat ) 
mod8.1 <- glm(prop18surv20 ~ region + winteravgmn, family = binomial(link = "logit") , data = dat ) 
mod9 <- glm(prop18surv20 ~ region + poly(springtotalsnow, 2), family = binomial(link = "logit") , data = dat ) 
mod9.1 <- glm(prop18surv20 ~ region + springtotalsnow, family = binomial(link = "logit") , data = dat ) 

aictab(list(canopycont2 = mod1, canopycont = mod1.1, f_b2 = mod2, f_b = mod2.1, whc2 = mod3, whc = mod3.1, 
            c_n2 = mod4, c_n = mod4.1, winteravgmn2 = mod8,
            winteravgmn = mod8.1, springtotalsnow2 = mod9, springtotalsnow = mod9.1))

# Fit global model (= full model)
mod <- glm(prop18surv20 ~ region + poly(canopycont, 2) + poly(f_b, 2) + poly(whc, 2) + c_n + 
             winteravgmn + springtotalsnow, family = binomial(link = "logit"), data = dat,
           na.action = "na.fail") 

# Model diagnostics I
vif(mod) 

# # CONCLUSION: removed no params

# Variable selection with dredge
dd <- dredge(mod, evaluate = T, rank = 'AICc', extra = 'R^2', m.lim = c(0, mm))

# Look at results
subset(dd, delta <= 2)

# Save output to use in microhab_suitab.R
saveRDS(subset(dd, delta <= 2), 'outputs/candmods_2yr_sorsit.rds') 

avgm <- model.avg(get.models(dd, subset = delta <= 2))
saveRDS(avgm, 'outputs/avgm_2yr_sorsit.rds') ###*** change here


# Model diagnostics II
mod <- glm(prop18surv20 ~ poly(canopycont, 2) + poly(f_b, 2) + poly(whc, 2) + region, 
           family = binomial(link = "logit"), data = dat) #fit best model

plot(fitted(mod),resid(mod)) #homogeneity of variance

lag.plot(resid(mod),diag=FALSE,do.lines=FALSE) #independence of residuals to neighboring residuals

deviance(mod)/df.residual(mod) #dispersion should = 1 for variance to equal the mean

coef(lm(with(dat,prop18surv20)~fitted(mod))) #link function appropriate if slope = 1

## CONCLUSION: 
# homogeneity of variance - met
# independence of residuals - met
# variance equals the mean - not met
# link function - not appropriate




####################################################################################################

# # TELGRA - 1 & 2 YR # # *** change here
# Family: Binomial

####################################################################################################

# Data - 1 YEAR SEEDLING SURVIVAL
dat <- census8.sdlg1yr %>%
  dplyr::filter(species == 'TELGRA') ###*** change here

# Identify how many parameters to use in model selection
mm <- round(nrow(dat)/10, 0) #rule of thumb: ~ 10 data points/parameter
mm

# Choose quadratic effects (<= 2 delta AICc)
mod1 <- glm(surv1yr ~ yr1 + region + poly(canopycont, 2), family = binomial(link = "logit"), data = dat ) #for AICc comparison
mod1.1 <- glm(surv1yr ~ yr1 + region + canopycont, family = binomial(link = "logit") , data = dat )
mod2 <- glm(surv1yr ~ yr1 + region + poly(f_b, 2), family = binomial(link = "logit") , data = dat ) 
mod2.1 <- glm(surv1yr ~ yr1 + region + f_b, family = binomial(link = "logit") , data = dat ) 
mod3 <- glm(surv1yr ~ yr1 + region + poly(whc, 2), family = binomial(link = "logit") , data = dat ) 
mod3.1 <- glm(surv1yr ~ yr1 + region + whc, family = binomial(link = "logit") , data = dat ) 
mod4 <- glm(surv1yr ~ yr1 + region + poly(c_n, 2), family = binomial(link = "logit") , data = dat ) 
mod4.1 <- glm(surv1yr ~ yr1 + region + c_n, family = binomial(link = "logit") , data = dat )
mod8 <- glm(surv1yr ~ yr1 + region + poly(winteravgmn, 2), family = binomial(link = "logit") , data = dat ) 
mod8.1 <- glm(surv1yr ~ yr1 + region + winteravgmn, family = binomial(link = "logit") , data = dat ) 
mod9 <- glm(surv1yr ~ yr1 + region + poly(springtotalsnow, 2), family = binomial(link = "logit") , data = dat ) 
mod9.1 <- glm(surv1yr ~ yr1 + region + springtotalsnow, family = binomial(link = "logit") , data = dat ) 

aictab(list(canopycont2 = mod1, canopycont = mod1.1, f_b2 = mod2, f_b = mod2.1, whc2 = mod3, whc = mod3.1, 
            c_n2 = mod4, c_n = mod4.1, winteravgmn2 = mod8,
            winteravgmn = mod8.1, springtotalsnow2 = mod9, springtotalsnow = mod9.1))

# Fit global model (= full model)
mod <- glm(surv1yr ~ yr1 + canopycont + f_b + whc + poly(c_n, 2) + 
             winteravgmn + springtotalsnow, family = binomial(link = "logit"), data = dat,
           na.action = "na.fail") 


# Model diagnostics I
vif(mod) 

# # CONCLUSION: removed region


# Variable selection with dredge
dd <- dredge(mod, evaluate = T, rank = 'AICc', extra = 'R^2', m.lim = c(0, mm))

# Look at results
subset(dd, delta <= 2)

# Save output to use in microhab_suitab.R
saveRDS(subset(dd, delta <= 2), 'outputs/candmods_1yr_telgra.rds') 

avgm <- model.avg(get.models(dd, subset = delta <= 2))
saveRDS(avgm, 'outputs/avgm_1yr_telgra.rds') ###*** change here


# Model diagnostics II
mod <- glm(surv1yr ~ canopycont + poly(c_n, 2) + springtotalsnow, 
           family = binomial(link = "logit"), data = dat) #fit best model

plot(fitted(mod),resid(mod)) #homogeneity of variance

lag.plot(resid(mod),diag=FALSE,do.lines=FALSE) #independence of residuals to neighboring residuals

deviance(mod)/df.residual(mod) #dispersion should = 1 for variance to equal the mean

coef(lm(with(dat,surv1yr)~fitted(mod))) #link function appropriate if slope = 1

## CONCLUSION: 
# homogeneity of variance - met
# independence of residuals - met
# variance equals the mean - not met
# link function - appropriate


# Data - 2 YEAR SEEDLING SURVIVAL
dat <- census8.sdlg2yr %>%
  dplyr::filter(species == 'TELGRA') ###*** change here

# Identify how many parameters to use in model selection
mm <- round(nrow(dat)/10, 0) #rule of thumb: ~ 10 data points/parameter
mm

# Choose quadratic effects (<= 2 delta AICc)
mod1 <- glm(prop18surv20 ~ region + poly(canopycont, 2), family = binomial(link = "logit"), data = dat ) #for AICc comparison
mod1.1 <- glm(prop18surv20 ~ region + canopycont, family = binomial(link = "logit") , data = dat )
mod2 <- glm(prop18surv20 ~ region + poly(f_b, 2), family = binomial(link = "logit") , data = dat ) 
mod2.1 <- glm(prop18surv20 ~ region + f_b, family = binomial(link = "logit") , data = dat ) 
mod3 <- glm(prop18surv20 ~ region + poly(whc, 2), family = binomial(link = "logit") , data = dat ) 
mod3.1 <- glm(prop18surv20 ~ region + whc, family = binomial(link = "logit") , data = dat ) 
mod4 <- glm(prop18surv20 ~ region + poly(c_n, 2), family = binomial(link = "logit") , data = dat ) 
mod4.1 <- glm(prop18surv20 ~ region + c_n, family = binomial(link = "logit") , data = dat )
mod8 <- glm(prop18surv20 ~ region + poly(winteravgmn, 2), family = binomial(link = "logit") , data = dat ) 
mod8.1 <- glm(prop18surv20 ~ region + winteravgmn, family = binomial(link = "logit") , data = dat ) 
mod9 <- glm(prop18surv20 ~ region + poly(springtotalsnow, 2), family = binomial(link = "logit") , data = dat ) 
mod9.1 <- glm(prop18surv20 ~ region + springtotalsnow, family = binomial(link = "logit") , data = dat ) 

aictab(list(canopycont2 = mod1, canopycont = mod1.1, f_b2 = mod2, f_b = mod2.1, whc2 = mod3, whc = mod3.1, 
            c_n2 = mod4, c_n = mod4.1, winteravgmn2 = mod8,
            winteravgmn = mod8.1, springtotalsnow2 = mod9, springtotalsnow = mod9.1))

# Fit global model (= full model)
mod <- glm(prop18surv20 ~ canopycont + f_b + whc + c_n + 
             winteravgmn, family = binomial(link = "logit"), data = dat,
           na.action = "na.fail") 

# Model diagnostics I
vif(mod) 

# # CONCLUSION: removed region, springtotalsnow

# Variable selection with dredge
dd <- dredge(mod, evaluate = T, rank = 'AICc', extra = 'R^2', m.lim = c(0, mm))

# Look at results
subset(dd, delta <= 2)

# Save output to use in microhab_suitab.R
saveRDS(subset(dd, delta <= 2), 'outputs/candmods_2yr_telgra.rds') 

avgm <- model.avg(get.models(dd, subset = delta <= 2))
saveRDS(avgm, 'outputs/avgm_2yr_telgra.rds') ###*** change here


# Model diagnostics II
mod <- glm(prop18surv20 ~ f_b, 
           family = binomial(link = "logit"), data = dat) #fit best model

plot(fitted(mod),resid(mod)) #homogeneity of variance

lag.plot(resid(mod),diag=FALSE,do.lines=FALSE) #independence of residuals to neighboring residuals

deviance(mod)/df.residual(mod) #dispersion should = 1 for variance to equal the mean

coef(lm(with(dat,prop18surv20)~fitted(mod))) #link function appropriate if slope = 1

## CONCLUSION: 
# homogeneity of variance - met
# independence of residuals - met
# variance equals the mean - moderately met
# link function - appropriate




####################################################################################################

# # TOLMEN - 1 & 2 YR # # *** change here
# Family: Binomial

####################################################################################################

# Data - 1 YEAR SEEDLING SURVIVAL
dat <- census8.sdlg1yr %>%
  dplyr::filter(species == 'TOLMEN') ###*** change here

# Identify how many parameters to use in model selection
mm <- round(nrow(dat)/10, 0) #rule of thumb: ~ 10 data points/parameter
mm

# Choose quadratic effects (<= 2 delta AICc)
mod1 <- glm(surv1yr ~ yr1 + region + poly(canopycont, 2), family = binomial(link = "logit"), data = dat ) #for AICc comparison
mod1.1 <- glm(surv1yr ~ yr1 + region + canopycont, family = binomial(link = "logit") , data = dat )
mod2 <- glm(surv1yr ~ yr1 + region + poly(f_b, 2), family = binomial(link = "logit") , data = dat ) 
mod2.1 <- glm(surv1yr ~ yr1 + region + f_b, family = binomial(link = "logit") , data = dat ) 
mod3 <- glm(surv1yr ~ yr1 + region + poly(whc, 2), family = binomial(link = "logit") , data = dat ) 
mod3.1 <- glm(surv1yr ~ yr1 + region + whc, family = binomial(link = "logit") , data = dat ) 
mod4 <- glm(surv1yr ~ yr1 + region + poly(c_n, 2), family = binomial(link = "logit") , data = dat ) 
mod4.1 <- glm(surv1yr ~ yr1 + region + c_n, family = binomial(link = "logit") , data = dat )
mod8 <- glm(surv1yr ~ yr1 + region + poly(winteravgmn, 2), family = binomial(link = "logit") , data = dat ) 
mod8.1 <- glm(surv1yr ~ yr1 + region + winteravgmn, family = binomial(link = "logit") , data = dat ) 
mod9 <- glm(surv1yr ~ yr1 + region + poly(springtotalsnow, 2), family = binomial(link = "logit") , data = dat ) 
mod9.1 <- glm(surv1yr ~ yr1 + region + springtotalsnow, family = binomial(link = "logit") , data = dat ) 

aictab(list(canopycont2 = mod1, canopycont = mod1.1, f_b2 = mod2, f_b = mod2.1, whc2 = mod3, whc = mod3.1, 
            c_n2 = mod4, c_n = mod4.1, winteravgmn2 = mod8,
            winteravgmn = mod8.1, springtotalsnow2 = mod9, springtotalsnow = mod9.1))

# Fit global model (= full model)
mod <- glm(surv1yr ~ yr1 + region + canopycont + f_b + poly(whc, 2) + c_n + 
             winteravgmn, family = binomial(link = "logit"), data = dat,
           na.action = "na.fail") 

# Model diagnostics I
vif(mod) 

# # CONCLUSION: removed poly springtotalsnow, springtotalsnow, poly winteravgmn


# Variable selection with dredge
dd <- dredge(mod, evaluate = T, rank = 'AICc', extra = 'R^2', m.lim = c(0, mm))

# Look at results
subset(dd, delta <= 2)

# Save output to use in microhab_suitab.R
saveRDS(subset(dd, delta <= 2), 'outputs/candmods_1yr_tolmen.rds') 

avgm <- model.avg(get.models(dd, subset = delta <= 2))
saveRDS(avgm, 'outputs/avgm_1yr_tolmen.rds') ###*** change here


# Model diagnostics II
mod <- glm(surv1yr ~ c_n + poly(whc, 2) + winteravgmn, 
           family = binomial(link = "logit"), data = dat) #fit best model

plot(fitted(mod),resid(mod)) #homogeneity of variance

lag.plot(resid(mod),diag=FALSE,do.lines=FALSE) #independence of residuals to neighboring residuals

deviance(mod)/df.residual(mod) #dispersion should = 1 for variance to equal the mean

coef(lm(with(dat,surv1yr)~fitted(mod))) #link function appropriate if slope = 1

## CONCLUSION: 
# homogeneity of variance - met
# independence of residuals - met
# variance equals the mean - not met
# link function - appropriate


# Data - 2 YEAR SEEDLING SURVIVAL
dat <- census8.sdlg2yr %>%
  dplyr::filter(species == 'TOLMEN') ###*** change here

# Identify how many parameters to use in model selection
mm <- round(nrow(dat)/10, 0) #rule of thumb: ~ 10 data points/parameter
mm

# Choose quadratic effects (<= 2 delta AICc)
mod1 <- glm(prop18surv20 ~ region + poly(canopycont, 2), family = binomial(link = "logit"), data = dat ) #for AICc comparison
mod1.1 <- glm(prop18surv20 ~ region + canopycont, family = binomial(link = "logit") , data = dat )
mod2 <- glm(prop18surv20 ~ region + poly(f_b, 2), family = binomial(link = "logit") , data = dat ) 
mod2.1 <- glm(prop18surv20 ~ region + f_b, family = binomial(link = "logit") , data = dat ) 
mod3 <- glm(prop18surv20 ~ region + poly(whc, 2), family = binomial(link = "logit") , data = dat ) 
mod3.1 <- glm(prop18surv20 ~ region + whc, family = binomial(link = "logit") , data = dat ) 
mod4 <- glm(prop18surv20 ~ region + poly(c_n, 2), family = binomial(link = "logit") , data = dat ) 
mod4.1 <- glm(prop18surv20 ~ region + c_n, family = binomial(link = "logit") , data = dat )
mod8 <- glm(prop18surv20 ~ region + poly(winteravgmn, 2), family = binomial(link = "logit") , data = dat ) 
mod8.1 <- glm(prop18surv20 ~ region + winteravgmn, family = binomial(link = "logit") , data = dat ) 
mod9 <- glm(prop18surv20 ~ region + poly(springtotalsnow, 2), family = binomial(link = "logit") , data = dat ) 
mod9.1 <- glm(prop18surv20 ~ region + springtotalsnow, family = binomial(link = "logit") , data = dat ) 

aictab(list(canopycont2 = mod1, canopycont = mod1.1, f_b2 = mod2, f_b = mod2.1, whc2 = mod3, whc = mod3.1, 
            c_n2 = mod4, c_n = mod4.1, winteravgmn2 = mod8,
            winteravgmn = mod8.1, springtotalsnow2 = mod9, springtotalsnow = mod9.1))

# Fit global model (= full model)
mod <- glm(prop18surv20 ~ canopycont + f_b + whc + c_n + 
             winteravgmn, family = binomial(link = "logit"), data = dat,
           na.action = "na.fail") 

# Model diagnostics I
vif(mod) 

# # CONCLUSION: removed region, springtotalsnow

# Variable selection with dredge
dd <- dredge(mod, evaluate = T, rank = 'AICc', extra = 'R^2', m.lim = c(0, mm))

# Look at results
subset(dd, delta <= 2)

# Model diagnostics II
mod <- glm(prop18surv20 ~ c_n + f_b, 
           family = binomial(link = "logit"), data = dat) #fit best model

# Save output to use in microhab_suitab.R
saveRDS(subset(dd, delta <= 2), 'outputs/candmods_2yr_tolmen.rds') 

saveRDS(mod, 'outputs/mod_2yr_tolmen.rds') ###*** change here

plot(fitted(mod),resid(mod)) #homogeneity of variance

lag.plot(resid(mod),diag=FALSE,do.lines=FALSE) #independence of residuals to neighboring residuals

deviance(mod)/df.residual(mod) #dispersion should = 1 for variance to equal the mean

coef(lm(with(dat,prop18surv20)~fitted(mod))) #link function appropriate if slope = 1

## CONCLUSION: 
# homogeneity of variance - met
# independence of residuals - met
# variance equals the mean - not met
# link function - appropriate




####################################################################################################

# # LUPLAT - 1 & 2 YR # # *** change here
# Family: Binomial

####################################################################################################

# Data - 1 YEAR SEEDLING SURVIVAL
dat <- census8.sdlg1yr %>%
  dplyr::filter(species == 'LUPLAT') ###*** change here

# Identify how many parameters to use in model selection
mm <- round(nrow(dat)/10, 0) #rule of thumb: ~ 10 data points/parameter
mm

# Choose quadratic effects (<= 2 delta AICc)
mod1 <- glm(surv1yr ~ yr1 + region + poly(canopycont, 2), family = binomial(link = "logit"), data = dat ) #for AICc comparison
mod1.1 <- glm(surv1yr ~ yr1 + region + canopycont, family = binomial(link = "logit") , data = dat )
mod2 <- glm(surv1yr ~ yr1 + region + poly(f_b, 2), family = binomial(link = "logit") , data = dat ) 
mod2.1 <- glm(surv1yr ~ yr1 + region + f_b, family = binomial(link = "logit") , data = dat ) 
mod3 <- glm(surv1yr ~ yr1 + region + poly(whc, 2), family = binomial(link = "logit") , data = dat ) 
mod3.1 <- glm(surv1yr ~ yr1 + region + whc, family = binomial(link = "logit") , data = dat ) 
mod4 <- glm(surv1yr ~ yr1 + region + poly(c_n, 2), family = binomial(link = "logit") , data = dat ) 
mod4.1 <- glm(surv1yr ~ yr1 + region + c_n, family = binomial(link = "logit") , data = dat )
mod8 <- glm(surv1yr ~ yr1 + region + poly(winteravgmn, 2), family = binomial(link = "logit") , data = dat ) 
mod8.1 <- glm(surv1yr ~ yr1 + region + winteravgmn, family = binomial(link = "logit") , data = dat ) 
mod9 <- glm(surv1yr ~ yr1 + region + poly(springtotalsnow, 2), family = binomial(link = "logit") , data = dat ) 
mod9.1 <- glm(surv1yr ~ yr1 + region + springtotalsnow, family = binomial(link = "logit") , data = dat ) 

aictab(list(canopycont2 = mod1, canopycont = mod1.1, f_b2 = mod2, f_b = mod2.1, whc2 = mod3, whc = mod3.1, 
            c_n2 = mod4, c_n = mod4.1, winteravgmn2 = mod8,
            winteravgmn = mod8.1, springtotalsnow2 = mod9, springtotalsnow = mod9.1))

# Fit global model (= full model)
mod <- glm(surv1yr ~ yr1 + region + poly(canopycont, 2) + f_b + whc + c_n + 
             winteravgmn + springtotalsnow, family = binomial(link = "logit"), data = dat,
           na.action = "na.fail") 


# Model diagnostics I
vif(mod) 

# # CONCLUSION: removed no params


# Variable selection with dredge
dd <- dredge(mod, evaluate = T, rank = 'AICc', extra = 'R^2', m.lim = c(0, mm))

# Look at results
subset(dd, delta <= 2)

# Save output to use in microhab_suitab.R
saveRDS(subset(dd, delta <= 2), 'outputs/candmods_1yr_luplat.rds') 

avgm <- model.avg(get.models(dd, subset = delta <= 2))
saveRDS(avgm, 'outputs/avgm_1yr_luplat.rds') ###*** change here


# Model diagnostics II
mod <- glm(surv1yr ~ poly(canopycont, 2) + whc, 
           family = binomial(link = "logit"), data = dat) #fit best model

plot(fitted(mod),resid(mod)) #homogeneity of variance

lag.plot(resid(mod),diag=FALSE,do.lines=FALSE) #independence of residuals to neighboring residuals

deviance(mod)/df.residual(mod) #dispersion should = 1 for variance to equal the mean

coef(lm(with(dat,surv1yr)~fitted(mod))) #link function appropriate if slope = 1

## CONCLUSION: 
# homogeneity of variance - met
# independence of residuals - met
# variance equals the mean - not met
# link function - appropriate


# Data - 2 YEAR SEEDLING SURVIVAL
dat <- census8.sdlg2yr %>%
  dplyr::filter(species == 'LUPLAT') ###*** change here

# Identify how many parameters to use in model selection
mm <- round(nrow(dat)/10, 0) #rule of thumb: ~ 10 data points/parameter
mm

# Choose quadratic effects (<= 2 delta AICc)
mod1 <- glm(prop18surv20 ~ region + poly(canopycont, 2), family = binomial(link = "logit"), data = dat ) #for AICc comparison
mod1.1 <- glm(prop18surv20 ~ region + canopycont, family = binomial(link = "logit") , data = dat )
mod2 <- glm(prop18surv20 ~ region + poly(f_b, 2), family = binomial(link = "logit") , data = dat ) 
mod2.1 <- glm(prop18surv20 ~ region + f_b, family = binomial(link = "logit") , data = dat ) 
mod3 <- glm(prop18surv20 ~ region + poly(whc, 2), family = binomial(link = "logit") , data = dat ) 
mod3.1 <- glm(prop18surv20 ~ region + whc, family = binomial(link = "logit") , data = dat ) 
mod4 <- glm(prop18surv20 ~ region + poly(c_n, 2), family = binomial(link = "logit") , data = dat ) 
mod4.1 <- glm(prop18surv20 ~ region + c_n, family = binomial(link = "logit") , data = dat )
mod8 <- glm(prop18surv20 ~ region + poly(winteravgmn, 2), family = binomial(link = "logit") , data = dat ) 
mod8.1 <- glm(prop18surv20 ~ region + winteravgmn, family = binomial(link = "logit") , data = dat ) 
mod9 <- glm(prop18surv20 ~ region + poly(springtotalsnow, 2), family = binomial(link = "logit") , data = dat ) 
mod9.1 <- glm(prop18surv20 ~ region + springtotalsnow, family = binomial(link = "logit") , data = dat ) 

aictab(list(canopycont2 = mod1, canopycont = mod1.1, f_b2 = mod2, f_b = mod2.1, whc2 = mod3, whc = mod3.1, 
            c_n2 = mod4, c_n = mod4.1, winteravgmn2 = mod8,
            winteravgmn = mod8.1, springtotalsnow2 = mod9, springtotalsnow = mod9.1))

# Fit global model (= full model)
mod <- glm(prop18surv20 ~ region + canopycont + f_b + whc + c_n + 
             winteravgmn, family = binomial(link = "logit"), data = dat,
           na.action = "na.fail") 

# Model diagnostics I
vif(mod) 

# # CONCLUSION: removed poly(c_n), springtotalsnow

# Variable selection with dredge
dd <- dredge(mod, evaluate = T, rank = 'AICc', extra = 'R^2', m.lim = c(0, mm))

# Look at results
subset(dd, delta <= 2)

# Save output to use in microhab_suitab.R
saveRDS(subset(dd, delta <= 2), 'outputs/candmods_2yr_luplat.rds') 

avgm <- model.avg(get.models(dd, subset = delta <= 2))
saveRDS(avgm, 'outputs/avgm_2yr_luplat.rds') ###*** change here


# Model diagnostics II
mod <- glm(prop18surv20 ~ f_b, 
           family = binomial(link = "logit"), data = dat) #fit best model

plot(fitted(mod),resid(mod)) #homogeneity of variance

lag.plot(resid(mod),diag=FALSE,do.lines=FALSE) #independence of residuals to neighboring residuals

deviance(mod)/df.residual(mod) #dispersion should = 1 for variance to equal the mean

coef(lm(with(dat,prop18surv20)~fitted(mod))) #link function appropriate if slope = 1

## CONCLUSION: 
# homogeneity of variance - met
# independence of residuals - met
# variance equals the mean - not met
# link function - appropriate




####################################################################################################

# # ERIPER - 1 & 2 YR # # *** change here
# Family: Binomial

####################################################################################################

# Data - 1 YEAR SEEDLING SURVIVAL
dat <- census8.sdlg1yr %>%
  dplyr::filter(species == 'ERIPER') ###*** change here

# Identify how many parameters to use in model selection
mm <- round(nrow(dat)/10, 0) #rule of thumb: ~ 10 data points/parameter
mm

# Choose quadratic effects (<= 2 delta AICc)
mod1 <- glm(surv1yr ~ yr1 + region + poly(canopycont, 2), family = binomial(link = "logit"), data = dat ) #for AICc comparison
mod1.1 <- glm(surv1yr ~ yr1 + region + canopycont, family = binomial(link = "logit") , data = dat )
mod2 <- glm(surv1yr ~ yr1 + region + poly(f_b, 2), family = binomial(link = "logit") , data = dat ) 
mod2.1 <- glm(surv1yr ~ yr1 + region + f_b, family = binomial(link = "logit") , data = dat ) 
mod3 <- glm(surv1yr ~ yr1 + region + poly(whc, 2), family = binomial(link = "logit") , data = dat ) 
mod3.1 <- glm(surv1yr ~ yr1 + region + whc, family = binomial(link = "logit") , data = dat ) 
mod4 <- glm(surv1yr ~ yr1 + region + poly(c_n, 2), family = binomial(link = "logit") , data = dat ) 
mod4.1 <- glm(surv1yr ~ yr1 + region + c_n, family = binomial(link = "logit") , data = dat )
mod8 <- glm(surv1yr ~ yr1 + region + poly(winteravgmn, 2), family = binomial(link = "logit") , data = dat ) 
mod8.1 <- glm(surv1yr ~ yr1 + region + winteravgmn, family = binomial(link = "logit") , data = dat ) 
mod9 <- glm(surv1yr ~ yr1 + region + poly(springtotalsnow, 2), family = binomial(link = "logit") , data = dat ) 
mod9.1 <- glm(surv1yr ~ yr1 + region + springtotalsnow, family = binomial(link = "logit") , data = dat ) 

aictab(list(canopycont2 = mod1, canopycont = mod1.1, f_b2 = mod2, f_b = mod2.1, whc2 = mod3, whc = mod3.1, 
            c_n2 = mod4, c_n = mod4.1, winteravgmn2 = mod8,
            winteravgmn = mod8.1, springtotalsnow2 = mod9, springtotalsnow = mod9.1))

# Fit global model (= full model)
mod <- glm(surv1yr ~ yr1 + region + canopycont + f_b + whc + c_n + 
             winteravgmn + springtotalsnow, family = binomial(link = "logit"), data = dat,
           na.action = "na.fail") 

# Model diagnostics I
vif(mod) 

# # CONCLUSION: removed poly(winteravgmn)


# Variable selection with dredge
dd <- dredge(mod, evaluate = T, rank = 'AICc', extra = 'R^2', m.lim = c(0, mm))

# Look at results
subset(dd, delta <= 2)

# Save output to use in microhab_suitab.R
saveRDS(subset(dd, delta <= 2), 'outputs/candmods_1yr_eriper.rds') 

avgm <- model.avg(get.models(dd, subset = delta <= 2))
saveRDS(avgm, 'outputs/avgm_1yr_eriper.rds') ###*** change here


# Model diagnostics II
mod <- glm(surv1yr ~ yr1, 
           family = binomial(link = "logit"), data = dat) #fit best model

plot(fitted(mod),resid(mod)) #homogeneity of variance

lag.plot(resid(mod),diag=FALSE,do.lines=FALSE) #independence of residuals to neighboring residuals

deviance(mod)/df.residual(mod) #dispersion should = 1 for variance to equal the mean

coef(lm(with(dat,surv1yr)~fitted(mod))) #link function appropriate if slope = 1

## CONCLUSION: 
# homogeneity of variance - not met but not expected
# independence of residuals - met
# variance equals the mean - met
# link function - appropriate


# Data - 2 YEAR SEEDLING SURVIVAL
dat <- census8.sdlg2yr %>%
  dplyr::filter(species == 'ERIPER') ###*** change here

# Identify how many parameters to use in model selection
mm <- round(nrow(dat)/10, 0) #rule of thumb: ~ 10 data points/parameter
mm

# Choose quadratic effects (<= 2 delta AICc)
mod1 <- glm(prop18surv20 ~ region + poly(canopycont, 2), family = binomial(link = "logit"), data = dat ) #for AICc comparison
mod1.1 <- glm(prop18surv20 ~ region + canopycont, family = binomial(link = "logit") , data = dat )
mod2 <- glm(prop18surv20 ~ region + poly(f_b, 2), family = binomial(link = "logit") , data = dat ) 
mod2.1 <- glm(prop18surv20 ~ region + f_b, family = binomial(link = "logit") , data = dat ) 
mod3 <- glm(prop18surv20 ~ region + poly(whc, 2), family = binomial(link = "logit") , data = dat ) 
mod3.1 <- glm(prop18surv20 ~ region + whc, family = binomial(link = "logit") , data = dat ) 
mod4 <- glm(prop18surv20 ~ region + poly(c_n, 2), family = binomial(link = "logit") , data = dat ) 
mod4.1 <- glm(prop18surv20 ~ region + c_n, family = binomial(link = "logit") , data = dat )
mod8 <- glm(prop18surv20 ~ region + poly(winteravgmn, 2), family = binomial(link = "logit") , data = dat ) 
mod8.1 <- glm(prop18surv20 ~ region + winteravgmn, family = binomial(link = "logit") , data = dat ) 
mod9 <- glm(prop18surv20 ~ region + poly(springtotalsnow, 2), family = binomial(link = "logit") , data = dat ) 
mod9.1 <- glm(prop18surv20 ~ region + springtotalsnow, family = binomial(link = "logit") , data = dat ) 

aictab(list(canopycont2 = mod1, canopycont = mod1.1, f_b2 = mod2, f_b = mod2.1, whc2 = mod3, whc = mod3.1, 
            c_n2 = mod4, c_n = mod4.1, winteravgmn2 = mod8,
            winteravgmn = mod8.1, springtotalsnow2 = mod9, springtotalsnow = mod9.1))

# Fit global model (= full model)
mod <- glm(prop18surv20 ~ region + canopycont + f_b + c_n + 
             winteravgmn + springtotalsnow, family = binomial(link = "logit"), data = dat,
           na.action = "na.fail") 

# Model diagnostics I
vif(mod) 

# # CONCLUSION: removed poly(canopycont), whc

# Variable selection with dredge
dd <- dredge(mod, evaluate = T, rank = 'AICc', extra = 'R^2', m.lim = c(0, mm))

# Look at results
subset(dd, delta <= 2)

# Save output to use in microhab_suitab.R
saveRDS(subset(dd, delta <= 2), 'outputs/candmods_2yr_eriper.rds') 

avgm <- model.avg(get.models(dd, subset = delta <= 2))
saveRDS(avgm, 'outputs/avgm_2yr_eriper.rds') ###*** change here

# Model diagnostics II
mod <- glm(prop18surv20 ~ canopycont + region + winteravgmn, 
           family = binomial(link = "logit"), data = dat) #fit best model

plot(fitted(mod),resid(mod)) #homogeneity of variance

lag.plot(resid(mod),diag=FALSE,do.lines=FALSE) #independence of residuals to neighboring residuals

deviance(mod)/df.residual(mod) #dispersion should = 1 for variance to equal the mean

coef(lm(with(dat,prop18surv20)~fitted(mod))) #link function appropriate if slope = 1

## CONCLUSION: 
# homogeneity of variance - met
# independence of residuals - met
# variance equals the mean - not met
# link function - appropriate




####################################################################################################

# # VACPAR - 1 YR # # *** change here
# Family: Binomial

####################################################################################################

# Data - 1 YEAR SEEDLING SURVIVAL
dat <- census8.sdlg1yr %>%
  dplyr::filter(species == 'VACPAR') ###*** change here

# Identify how many parameters to use in model selection
mm <- round(nrow(dat)/10, 0) #rule of thumb: ~ 10 data points/parameter
mm

# Choose quadratic effects (<= 2 delta AICc)
mod1 <- glm(surv1yr ~ yr1 + region + poly(canopycont, 2), family = binomial(link = "logit"), data = dat ) #for AICc comparison
mod1.1 <- glm(surv1yr ~ yr1 + region + canopycont, family = binomial(link = "logit") , data = dat )
mod2 <- glm(surv1yr ~ yr1 + region + poly(f_b, 2), family = binomial(link = "logit") , data = dat ) 
mod2.1 <- glm(surv1yr ~ yr1 + region + f_b, family = binomial(link = "logit") , data = dat ) 
mod3 <- glm(surv1yr ~ yr1 + region + poly(whc, 2), family = binomial(link = "logit") , data = dat ) 
mod3.1 <- glm(surv1yr ~ yr1 + region + whc, family = binomial(link = "logit") , data = dat ) 
mod4 <- glm(surv1yr ~ yr1 + region + poly(c_n, 2), family = binomial(link = "logit") , data = dat ) 
mod4.1 <- glm(surv1yr ~ yr1 + region + c_n, family = binomial(link = "logit") , data = dat )
mod8 <- glm(surv1yr ~ yr1 + region + poly(winteravgmn, 2), family = binomial(link = "logit") , data = dat ) 
mod8.1 <- glm(surv1yr ~ yr1 + region + winteravgmn, family = binomial(link = "logit") , data = dat ) 
mod9 <- glm(surv1yr ~ yr1 + region + poly(springtotalsnow, 2), family = binomial(link = "logit") , data = dat ) 
mod9.1 <- glm(surv1yr ~ yr1 + region + springtotalsnow, family = binomial(link = "logit") , data = dat ) 

aictab(list(canopycont2 = mod1, canopycont = mod1.1, f_b2 = mod2, f_b = mod2.1, whc2 = mod3, whc = mod3.1, 
            c_n2 = mod4, c_n = mod4.1, winteravgmn2 = mod8,
            winteravgmn = mod8.1, springtotalsnow2 = mod9, springtotalsnow = mod9.1))

# Fit global model (= full model)
mod <- glm(surv1yr ~ yr1 + region + canopycont + f_b + whc + c_n + 
             winteravgmn + springtotalsnow, family = binomial(link = "logit"), data = dat,
           na.action = "na.fail") 

# Model diagnostics I
vif(mod) 

# # CONCLUSION: removed no params


# Variable selection with dredge
dd <- dredge(mod, evaluate = T, rank = 'AICc', extra = 'R^2', m.lim = c(0, mm))

# Look at results
subset(dd, delta <= 2)

# Save output to use in microhab_suitab.R
saveRDS(subset(dd, delta <= 2), 'outputs/candmods_1yr_vacpar.rds') 

avgm <- model.avg(get.models(dd, subset = delta <= 2))
saveRDS(avgm, 'outputs/avgm_1yr_vacpar.rds') ###*** change here


# Model diagnostics II
mod <- glm(surv1yr ~ c_n + yr1, 
           family = binomial(link = "logit"), data = dat) #fit best model

plot(fitted(mod),resid(mod)) #homogeneity of variance

lag.plot(resid(mod),diag=FALSE,do.lines=FALSE) #independence of residuals to neighboring residuals

deviance(mod)/df.residual(mod) #dispersion should = 1 for variance to equal the mean

coef(lm(with(dat,surv1yr)~fitted(mod))) #link function appropriate if slope = 1

## CONCLUSION: 
# homogeneity of variance - met
# independence of residuals - met
# variance equals the mean - not met
# link function - appropriate




####################################################################################################

# # VACDEL - 1 YR # # *** change here
# Family: Binomial

####################################################################################################

# Data - 1 YEAR SEEDLING SURVIVAL
dat <- census8.sdlg1yr %>%
  dplyr::filter(species == 'VACDEL') ###*** change here

# Identify how many parameters to use in model selection
mm <- round(nrow(dat)/10, 0) #rule of thumb: ~ 10 data points/parameter
mm

# Choose quadratic effects (<= 2 delta AICc)
mod1 <- glm(surv1yr ~ yr1 + region + poly(canopycont, 2), family = binomial(link = "logit"), data = dat ) #for AICc comparison
mod1.1 <- glm(surv1yr ~ yr1 + region + canopycont, family = binomial(link = "logit") , data = dat )
mod2 <- glm(surv1yr ~ yr1 + region + poly(f_b, 2), family = binomial(link = "logit") , data = dat ) 
mod2.1 <- glm(surv1yr ~ yr1 + region + f_b, family = binomial(link = "logit") , data = dat ) 
mod3 <- glm(surv1yr ~ yr1 + region + poly(whc, 2), family = binomial(link = "logit") , data = dat ) 
mod3.1 <- glm(surv1yr ~ yr1 + region + whc, family = binomial(link = "logit") , data = dat ) 
mod4 <- glm(surv1yr ~ yr1 + region + poly(c_n, 2), family = binomial(link = "logit") , data = dat ) 
mod4.1 <- glm(surv1yr ~ yr1 + region + c_n, family = binomial(link = "logit") , data = dat )
mod8 <- glm(surv1yr ~ yr1 + region + poly(winteravgmn, 2), family = binomial(link = "logit") , data = dat ) 
mod8.1 <- glm(surv1yr ~ yr1 + region + winteravgmn, family = binomial(link = "logit") , data = dat ) 
mod9 <- glm(surv1yr ~ yr1 + region + poly(springtotalsnow, 2), family = binomial(link = "logit") , data = dat ) 
mod9.1 <- glm(surv1yr ~ yr1 + region + springtotalsnow, family = binomial(link = "logit") , data = dat ) 

aictab(list(canopycont2 = mod1, canopycont = mod1.1, f_b2 = mod2, f_b = mod2.1, whc2 = mod3, whc = mod3.1, 
            c_n2 = mod4, c_n = mod4.1, winteravgmn2 = mod8,
            winteravgmn = mod8.1, springtotalsnow2 = mod9, springtotalsnow = mod9.1))

# Fit global model (= full model)
mod <- glm(surv1yr ~ yr1 + region + canopycont + f_b + whc + c_n + 
             springtotalsnow, family = binomial(link = "logit"), data = dat,
           na.action = "na.fail") 

# Model diagnostics I
vif(mod) 

# # CONCLUSION: removed winteravgmn, poly(c_n)


# Variable selection with dredge
dd <- dredge(mod, evaluate = T, rank = 'AICc', extra = 'R^2', m.lim = c(0, mm))

# Look at results
subset(dd, delta <= 2)

# Save output to use in microhab_suitab.R
saveRDS(subset(dd, delta <= 2), 'outputs/candmods_1yr_vacdel.rds') 

avgm <- model.avg(get.models(dd, subset = delta <= 2))
saveRDS(avgm, 'outputs/avgm_1yr_vacdel.rds') ###*** change here


# Model diagnostics II
mod <- glm(surv1yr ~ f_b, 
           family = binomial(link = "logit"), data = dat) #fit best model

plot(fitted(mod),resid(mod)) #homogeneity of variance

lag.plot(resid(mod),diag=FALSE,do.lines=FALSE) #independence of residuals to neighboring residuals

deviance(mod)/df.residual(mod) #dispersion should = 1 for variance to equal the mean

coef(lm(with(dat,surv1yr)~fitted(mod))) #link function appropriate if slope = 1

## CONCLUSION: 
# homogeneity of variance - met
# independence of residuals - met
# variance equals the mean - not met
# link function - appropriate



