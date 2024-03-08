# Aims:
# 1. Fit models for rel_rec > 0 with negative binomial GLM
# 3. Model averaging with MuMIn for species-level models (continuous recruitment data)

# Date created: 15 June 2023
# Date updated: 15 June 2023

# # LIBRARIES # #
library(tidyverse)
library(dplyr)
library(MuMIn)
library(AICcmodavg)
library(DHARMa)
library(car)

rm(list=ls()) 


# # INPUT FILES # #
load('data/census_spp.RData') #adds relative recruitment to census.tms (microhab_analysis_prep.R)


# # OUTPUT FILES # #

# Continuous Recruitment Model averaging and candidate models
avgm <- read_rds('outputs/avgm_cont_mahner.rds')
candmods <- read_rds('outputs/candmods_cont_mahner.rds')

avgm <- read_rds('outputs/mod_cont_erilan.rds') #best model
candmods <- read_rds('outputs/candmods_cont_erilan.rds')

avgm <- read_rds('outputs/avgm_cont_vacpar.rds')
candmods <- read_rds('outputs/candmods_cont_vacpar.rds')

avgm <- read_rds('outputs/avgm_cont_luplat.rds')
candmods <- read_rds('outputs/candmods_cont_luplat.rds')

avgm <- read_rds('outputs/avgm_cont_abilas.rds') 
candmods <- read_rds('outputs/candmods_cont_abilas.rds')

avgm <- read_rds('outputs/avgm_cont_aneocc.rds')
candmods <- read_rds('outputs/candmods_cont_aneocc.rds')

avgm <- read_rds('outputs/avgm_cont_eriper.rds')
candmods <- read_rds('outputs/candmods_cont_eriper.rds')

avgm <- read_rds('outputs/avgm_cont_piceng.rds') 
candmods <- read_rds('outputs/candmods_cont_piceng.rds')

avgm <- read_rds('outputs/avgm_cont_ruburs.rds')
candmods <- read_rds('outputs/candmods_cont_ruburs.rds')

avgm <- read_rds('outputs/avgm_cont_sorsit.rds')
candmods <- read_rds('outputs/candmods_cont_sorsit.rds')

avgm <- read_rds('outputs/avgm_cont_telgra.rds') 
candmods <- read_rds('outputs/candmods_cont_telgra.rds')

avgm <- read_rds('outputs/avgm_cont_tolmen.rds')
candmods <- read_rds('outputs/candmods_cont_tolmen.rds')

avgm <- read_rds('outputs/avgm_cont_vacdel.rds')
candmods <- read_rds('outputs/candmods_cont_vacdel.rds')

avgm <- read_rds('outputs/avgm_cont_mahaqu.rds')
candmods <- read_rds('outputs/candmods_cont_mahaqu.rds')




####################################################################################################

# # SET UP DATA FOR MODELS (no TOMST data) # # 

####################################################################################################

# Data
load('data/census_spp.RData') #adds relative recruitment to census.tms (microhab_analysis_prep.R)

# Species with indiv in >= 8 plots and rel_rec > 0
census8 <- census.spp %>% 
  group_by(species) %>% 
  filter(species == "LUPLAT" | species == "ABILAS" | species == "ANEOCC" |
           species == "ERIPER" | species == "PICENG" |
           species == "RUBURS" | species == "SORSIT" | species == "TELGRA" |
           species == "TOLMEN" | species == "VACDEL" | species == "VACPAR" |
           species == "MAHNER" | species == "ERILAN" | species == "MAHAQU") %>% 
  filter(rel_rec > 0)

# Scale variables (without TOMST)
census8 <- census8 %>%
  
  filter_at(vars(rel_rec, region, canopycont, f_b, whc, c_n, winteravgmn, springtotalsnow),
            all_vars(!is.na(.))) %>% #remove NAs
  
  mutate_at(vars(canopycont, f_b, whc, c_n, winteravgmn, springtotalsnow), ~scale(.x)[, 1])


# # Full model (with TOMST)
# mod <- lm(rel_rec ~ region + canopycont + f_b + whc + c_n + t1_avgmax + tmoist_avgmin + 
#             winteravgmn + springtotalsnow, data = dat) 
# 
# 
# # Full model (without TOMST)
# mod <- lm(rel_rec ~ region + canopycont + f_b + whc + c_n + winteravgmn + springtotalsnow, data = dat) 




####################################################################################################

# # CHECK GERMINATION AT SITES WITH MISSING MICROHAB DATA # # 

####################################################################################################

# Data
load('data/census_spp.RData') #adds relative recruitment to census.tms (microhab_analysis_prep.R)

census8 <- census.spp %>% 
  group_by(species) %>% 
  filter(species == "LUPLAT" | species == "ABILAS" | species == "ANEOCC" |
           species == "ERIPER" | species == "PICENG" |
           species == "RUBURS" | species == "SORSIT" | species == "TELGRA" |
           species == "TOLMEN" | species == "VACDEL" | species == "VACPAR" |
           species == "MAHNER" | species == "ERILAN" | species == "MAHAQU")

# Which variables have many NAs?
census8 %>% 
  select(canopycont, f_b, whc, c_n, tmoist_avgmin, t1_avgmax, t3_avgmin, winteravgmn, springtotalsnow) %>% 
  summary()

# # Conclusion: T3 > T1, Tmoist

miss.base <- census8 %>% 
  select(new_all, species) %>% #baseline of all recruitment data
  group_by(species) %>% 
  summarise(min = min(new_all), mean(new_all), max = max(new_all))

miss.t3 <- census8 %>% 
  filter(is.na(t3_avgmin)) %>% #NA values for TOMST T3 sensor
  select(new_all, species) %>% 
  group_by(species) %>% 
  summarise(min = min(new_all), mean(new_all), max = max(new_all))

# t3_avgmin NA values don't affect: ERILAN, MAHNER

miss.t1 <- census8 %>% 
  filter(is.na(t1_avgmax)) %>% #NA values for TOMST T1 sensor (NA for T1 = NA for Tmoist)
  select(new_all, species) %>% 
  group_by(species) %>% 
  summarise(min = min(new_all), mean(new_all), max = max(new_all))

# t1_avgmax NA values don't affect: ERILAN, MAHNER 

# # CONCLUSION: run only ERILAN & MAHNER with TOMST variables and remove from all others (3.5.23)




####################################################################################################

# # MAHNER - MODEL AVERAGING (with TOMST data) # #
# Family: Binomial

####################################################################################################

# Resources: help page for 'dredge', recruitment_optima_analysis.R (KG) for model diagnostics

# Data
dat <- census.spp %>% 
  
  filter_at(vars(rel_rec, region, canopycont, f_b, whc, c_n, tmoist_avgmin,
                 t1_avgmax, winteravgmn, springtotalsnow), all_vars(!is.na(.))) %>% #remove NAs
  
  mutate_at(vars(canopycont, f_b, whc, c_n, tmoist_avgmin,
                 t1_avgmax, winteravgmn, springtotalsnow), ~scale(.x)[, 1]) %>% #scale variables
  
  dplyr::filter(species == 'MAHNER') %>% 
  dplyr::filter(rel_rec > 0)

# Identify how many parameters to use in model selection
mm <- round(nrow(dat)/10, 0) #rule of thumb: ~ 10 data points/parameter
mm

# Choose quadratic effects (<= 2 delta AICc)
mod1 <- glm(rel_rec ~ region + poly(canopycont, 2), family = binomial(link = "logit"), data = dat )
mod1.1 <- glm(rel_rec ~ region + canopycont, family = binomial(link = "logit") , data = dat )
mod2 <- glm(rel_rec ~ region + poly(f_b, 2), family = binomial(link = "logit") , data = dat ) 
mod2.1 <- glm(rel_rec ~ region + f_b, family = binomial(link = "logit") , data = dat ) 
mod3 <- glm(rel_rec ~ region + poly(whc, 2), family = binomial(link = "logit") , data = dat ) 
mod3.1 <- glm(rel_rec ~ region + whc, family = binomial(link = "logit") , data = dat ) 
mod4 <- glm(rel_rec ~ region + poly(c_n, 2), family = binomial(link = "logit") , data = dat ) 
mod4.1 <- glm(rel_rec ~ region + c_n, family = binomial(link = "logit") , data = dat )
mod5 <- glm(rel_rec ~ region + poly(t1_avgmax, 2), family = binomial(link = "logit") , data = dat ) 
mod5.1 <- glm(rel_rec ~ region + t1_avgmax, family = binomial(link = "logit") , data = dat ) 
mod7 <- glm(rel_rec ~ region + poly(tmoist_avgmin, 2), family = binomial(link = "logit") , data = dat ) 
mod7.1 <- glm(rel_rec ~ region + tmoist_avgmin, family = binomial(link = "logit") , data = dat ) 
mod8 <- glm(rel_rec ~ region + poly(winteravgmn, 2), family = binomial(link = "logit") , data = dat ) 
mod8.1 <- glm(rel_rec ~ region + winteravgmn, family = binomial(link = "logit") , data = dat ) 
mod9 <- glm(rel_rec ~ region + poly(springtotalsnow, 2), family = binomial(link = "logit") , data = dat ) 
mod9.1 <- glm(rel_rec ~ region + springtotalsnow, family = binomial(link = "logit") , data = dat ) 

aictab(list(canopycont2 = mod1, canopycont = mod1.1, f_b2 = mod2, f_b = mod2.1, whc2 = mod3, whc = mod3.1, 
            c_n2 = mod4, c_n = mod4.1, t1_avgmax2 = mod5, t1_avgmax = mod5.1, 
            tmoist_avgmin2 = mod7, tmoist_avgmin = mod7.1, winteravgmn2 = mod8,
            winteravgmn = mod8.1, springtotalsnow2 = mod9, springtotalsnow = mod9.1))


# Fit global model (= full model)
mod <- glm(rel_rec ~ region + canopycont + f_b + whc + c_n + t1_avgmax + tmoist_avgmin + 
          winteravgmn + springtotalsnow, family = binomial(link = "logit"), data = dat,
           na.action = "na.fail") 

# Model diagnostics I
# vif(mod) #not needed for univariate model

# # CONCLUSION: no params taken out

# Variable selection with dredge
dd <- dredge(mod, evaluate = T, rank = 'AICc', extra = 'R^2', m.lim = c(0, mm))

# Save output to use in microhab_suitab.R
saveRDS(subset(dd, delta <= 2), 'outputs/candmods_cont_mahner.rds') 

avgm.cont.mahner <- model.avg(get.models(dd, subset = delta <= 2))
saveRDS(avgm.cont.mahner, 'outputs/avgm_cont_mahner.rds')

# Look at results
subset(dd, delta <= 2)


# Model diagnostics II
mod <- glm(rel_rec ~ whc, family = binomial(link = "logit") , data = dat ) #fit best model

plot(fitted(mod),resid(mod)) #homogeneity of variance

lag.plot(resid(mod),diag=FALSE,do.lines=FALSE) #independence of residuals to neighboring residuals

deviance(mod)/df.residual(mod) #dispersion should = 1 for variance to equal the mean

coef(lm(with(dat,rel_rec)~fitted(mod))) #link function appropriate if slope = 1

## CONCLUSION: 
# homogeneity of variance - met
# independence of residuals - met
# variance equals the mean - not met
# link function - appropriate



####################################################################################################

# # ERILAN - BEST MODEL (with TOMST data) # #
# Family: Binomial

####################################################################################################

# Data
dat <- census.spp %>% 
  
  filter_at(vars(rel_rec, region, canopycont, f_b, whc, c_n, tmoist_avgmin,
                 t1_avgmax, winteravgmn, springtotalsnow), all_vars(!is.na(.))) %>% #remove NAs
  
  mutate_at(vars(canopycont, f_b, whc, c_n, tmoist_avgmin,
                 t1_avgmax, winteravgmn, springtotalsnow), ~scale(.x)[, 1]) %>% #scale variables
  
  dplyr::filter(species == 'ERILAN') %>% 
  dplyr::filter(rel_rec > 0)

# Identify how many parameters to use in model selection
mm <- round(nrow(dat)/10, 0) #rule of thumb: ~ 10 data points/parameter
mm

# Choose quadratic effects (<= 2 delta AICc)
mod1 <- glm(rel_rec ~ poly(canopycont, 2), family = binomial(link = "logit"), data = dat )
mod1.1 <- glm(rel_rec ~ canopycont, family = binomial(link = "logit") , data = dat )
mod2 <- glm(rel_rec ~ poly(f_b, 2), family = binomial(link = "logit") , data = dat ) 
mod2.1 <- glm(rel_rec ~ f_b, family = binomial(link = "logit") , data = dat ) 
mod3 <- glm(rel_rec ~ poly(whc, 2), family = binomial(link = "logit") , data = dat ) 
mod3.1 <- glm(rel_rec ~ whc, family = binomial(link = "logit") , data = dat ) 
mod4 <- glm(rel_rec ~ poly(c_n, 2), family = binomial(link = "logit") , data = dat ) 
mod4.1 <- glm(rel_rec ~ c_n, family = binomial(link = "logit") , data = dat )
mod5 <- glm(rel_rec ~ poly(t1_avgmax, 2), family = binomial(link = "logit") , data = dat ) 
mod5.1 <- glm(rel_rec ~ t1_avgmax, family = binomial(link = "logit") , data = dat ) 
mod7 <- glm(rel_rec ~ poly(tmoist_avgmin, 2), family = binomial(link = "logit") , data = dat ) 
mod7.1 <- glm(rel_rec ~ tmoist_avgmin, family = binomial(link = "logit") , data = dat ) 
mod8 <- glm(rel_rec ~ poly(winteravgmn, 2), family = binomial(link = "logit") , data = dat ) 
mod8.1 <- glm(rel_rec ~ winteravgmn, family = binomial(link = "logit") , data = dat ) 
mod9 <- glm(rel_rec ~ poly(springtotalsnow, 2), family = binomial(link = "logit") , data = dat ) 
mod9.1 <- glm(rel_rec ~ springtotalsnow, family = binomial(link = "logit") , data = dat ) 

aictab(list(canopycont2 = mod1, canopycont = mod1.1, f_b2 = mod2, f_b = mod2.1, whc2 = mod3, whc = mod3.1, 
            c_n2 = mod4, c_n = mod4.1, t1_avgmax2 = mod5, t1_avgmax = mod5.1, 
            tmoist_avgmin2 = mod7, tmoist_avgmin = mod7.1, winteravgmn2 = mod8,
            winteravgmn = mod8.1, springtotalsnow2 = mod9, springtotalsnow = mod9.1))


# Fit global model (= full model)
mod <- glm(rel_rec ~ canopycont + f_b + whc + c_n + t1_avgmax + tmoist_avgmin + 
          winteravgmn + springtotalsnow, family = binomial(link = "logit"), data = dat,
           na.action = "na.fail") 

# Model diagnostics I
# vif(mod) #no need for univariate models

# # CONCLUSION: took out region (only 1 region germinated)


# Variable selection with dredge
dd <- dredge(mod, evaluate = T, rank = 'AICc', extra = 'R^2', m.lim = c(0, mm))

# Look at results
subset(dd, delta <= 2)


# Model diagnostics II
mod <- glm(formula = rel_rec ~ 1, family = binomial(link = "logit"), data = dat) #fit best model

# Save model to use in microhab_suitab.R
saveRDS(subset(dd, delta <= 2), 'outputs/candmods_cont_erilan.rds') 

saveRDS(mod, 'outputs/mod_cont_erilan.rds') ###*** change here


plot(fitted(mod),resid(mod)) #homogeneity of variance

lag.plot(resid(mod),diag=FALSE,do.lines=FALSE) #independence of residuals to neighboring residuals

deviance(mod)/df.residual(mod) #dispersion should = 1 for variance to equal the mean

coef(lm(with(dat,rel_rec)~fitted(mod))) #link function appropriate if slope = 1

## CONCLUSION: 
# homogeneity of variance - not met 
# independence of residuals - met
# variance equals the mean - not met
# link function - couldn't be tested with intercept-only model




####################################################################################################

# # VACPAR - MODEL AVERAGING # #
# Family: Binomial

####################################################################################################

# Data
dat <- census8 %>% 
  dplyr::filter(species == 'VACPAR') 

# Identify how many parameters to use in model selection
mm <- round(nrow(dat)/10, 0) #rule of thumb: ~ 10 data points/parameter
mm


# Choose quadratic effects (< 2 delta AICc)
mod1 <- glm(rel_rec ~ region + poly(canopycont, 2), family = binomial(link = "logit"), data = dat )
mod1.1 <- glm(rel_rec ~ region + canopycont, family = binomial(link = "logit") , data = dat )
mod2 <- glm(rel_rec ~ region + poly(f_b, 2), family = binomial(link = "logit") , data = dat ) 
mod2.1 <- glm(rel_rec ~ region + f_b, family = binomial(link = "logit") , data = dat ) 
mod3 <- glm(rel_rec ~ region + poly(whc, 2), family = binomial(link = "logit") , data = dat ) 
mod3.1 <- glm(rel_rec ~ region + whc, family = binomial(link = "logit") , data = dat ) 
mod4 <- glm(rel_rec ~ region + poly(c_n, 2), family = binomial(link = "logit") , data = dat ) 
mod4.1 <- glm(rel_rec ~ region + c_n, family = binomial(link = "logit") , data = dat )
mod8 <- glm(rel_rec ~ region + poly(winteravgmn, 2), family = binomial(link = "logit") , data = dat ) 
mod8.1 <- glm(rel_rec ~ region + winteravgmn, family = binomial(link = "logit") , data = dat ) 
mod9 <- glm(rel_rec ~ region + poly(springtotalsnow, 2), family = binomial(link = "logit") , data = dat ) 
mod9.1 <- glm(rel_rec ~ region + springtotalsnow, family = binomial(link = "logit") , data = dat ) 

aictab(list(canopycont2 = mod1, canopycont = mod1.1, f_b2 = mod2, f_b = mod2.1, whc2 = mod3, whc = mod3.1, 
            c_n2 = mod4, c_n = mod4.1, winteravgmn2 = mod8,
            winteravgmn = mod8.1, springtotalsnow2 = mod9, springtotalsnow = mod9.1))


# Fit global model (= full model)
mod <- glm(rel_rec ~ region + canopycont + f_b + whc + poly(c_n, 2) +
             winteravgmn + poly(springtotalsnow, 2), family = binomial(link = "logit"), data = dat,
           na.action = "na.fail") 

# Model diagnostics I
vif(mod) #check for collinearity > 5

# # CONCLUSION: removed no params


# Variable selection with dredge
dd <- dredge(mod, evaluate = T, rank = 'AICc', extra = 'R^2', m.lim = c(0, mm)) 

# Look at results with delta AICc < 2
subset(dd, delta <= 2)

# Save output to use in microhab_suitab.R
saveRDS(subset(dd, delta <= 2), 'outputs/candmods_cont_vacpar.rds') 

avgm <- model.avg(get.models(dd, subset = delta <= 2))
saveRDS(avgm, 'outputs/avgm_cont_vacpar.rds') ###*** change here


 # Model diagnostics II
mod <- glm(rel_rec ~ whc,
           family = binomial(link = "logit"), data = dat) #fit best model

plot(fitted(mod),resid(mod)) #homogeneity of variance

lag.plot(resid(mod),diag=FALSE,do.lines=FALSE) #independence of residuals to neighboring residuals

deviance(mod)/df.residual(mod) #dispersion should = 1 for variance to equal the mean

coef(lm(with(dat,rel_rec)~fitted(mod))) #link function appropriate if slope = 1

## CONCLUSION: 
# homogeneity of variance - met
# independence of residuals - met
# variance equals the mean - not met
# link function - moderately appropriate




####################################################################################################

# # LUPLAT - MODEL AVERAGING # #
# Family: Binomial

####################################################################################################

# Data
dat <- census8 %>% 
  dplyr::filter(species == 'LUPLAT')

# Identify how many parameters to use in model selection
mm <- round(nrow(dat)/10, 0) #rule of thumb: ~ 10 data points/parameter
mm

# Choose quadratic effects (< 2 delta AICc)
mod1 <- glm(rel_rec ~ region + poly(canopycont, 2), family = binomial(link = "logit"), data = dat )
mod1.1 <- glm(rel_rec ~ region + canopycont, family = binomial(link = "logit") , data = dat )
mod2 <- glm(rel_rec ~ region + poly(f_b, 2), family = binomial(link = "logit") , data = dat ) 
mod2.1 <- glm(rel_rec ~ region + f_b, family = binomial(link = "logit") , data = dat ) 
mod3 <- glm(rel_rec ~ region + poly(whc, 2), family = binomial(link = "logit") , data = dat ) 
mod3.1 <- glm(rel_rec ~ region + whc, family = binomial(link = "logit") , data = dat ) 
mod4 <- glm(rel_rec ~ region + poly(c_n, 2), family = binomial(link = "logit") , data = dat ) 
mod4.1 <- glm(rel_rec ~ region + c_n, family = binomial(link = "logit") , data = dat )
mod8 <- glm(rel_rec ~ region + poly(winteravgmn, 2), family = binomial(link = "logit") , data = dat ) 
mod8.1 <- glm(rel_rec ~ region + winteravgmn, family = binomial(link = "logit") , data = dat ) 
mod9 <- glm(rel_rec ~ region + poly(springtotalsnow, 2), family = binomial(link = "logit") , data = dat ) 
mod9.1 <- glm(rel_rec ~ region + springtotalsnow, family = binomial(link = "logit") , data = dat ) 

aictab(list(canopycont2 = mod1, canopycont = mod1.1, f_b2 = mod2, f_b = mod2.1, whc2 = mod3, whc = mod3.1, 
            c_n2 = mod4, c_n = mod4.1, winteravgmn2 = mod8,
            winteravgmn = mod8.1, springtotalsnow2 = mod9, springtotalsnow = mod9.1))


# Fit global model (= full model)
mod <- glm(rel_rec ~ region + canopycont + f_b + whc + c_n + winteravgmn + springtotalsnow, 
           family = binomial(link = "logit"), data = dat, na.action = "na.fail") 

# Model diagnostics
vif(mod) #check for collinearity > 5

# # CONCLUSION: removed no params

# Variable selection with dredge
dd <- dredge(mod, evaluate = T, rank = 'AICc', extra = 'R^2', m.lim = c(0, mm))

# Look at results
subset(dd, delta <= 2)

# Save output to use in microhab_suitab.R
saveRDS(subset(dd, delta <= 2), 'outputs/candmods_cont_luplat.rds') 

avgm <- model.avg(get.models(dd, subset = delta <= 2))
saveRDS(avgm, 'outputs/avgm_cont_luplat.rds') ###*** change here


# Model diagnostics II
mod <- glm(rel_rec ~ winteravgmn,
           family = binomial(link = "logit"), data = dat) #fit best model

plot(fitted(mod),resid(mod)) #homogeneity of variance

lag.plot(resid(mod),diag=FALSE,do.lines=FALSE) #independence of residuals to neighboring residuals

deviance(mod)/df.residual(mod) #dispersion should = 1 for variance to equal the mean

coef(lm(with(dat,rel_rec)~fitted(mod))) #link function appropriate if slope = 1

## CONCLUSION: 
# homogeneity of variance - not met
# independence of residuals - met
# variance equals the mean - not met
# link function - appropriate




####################################################################################################

# # ABILAS - MODEL AVERAGING # #
# Family: Binomial

####################################################################################################

# Data
dat <- census8 %>% 
  dplyr::filter(species == 'ABILAS') 

# Identify how many parameters to use in model selection
mm <- round(nrow(dat)/10, 0) #rule of thumb: ~ 10 data points/parameter
mm

# Choose quadratic effects (<= 2 delta AICc)
mod1 <- glm(rel_rec ~ region + poly(canopycont, 2), family = binomial(link = "logit"), data = dat )
mod1.1 <- glm(rel_rec ~ region + canopycont, family = binomial(link = "logit") , data = dat )
mod2 <- glm(rel_rec ~ region + poly(f_b, 2), family = binomial(link = "logit") , data = dat ) 
mod2.1 <- glm(rel_rec ~ region + f_b, family = binomial(link = "logit") , data = dat ) 
mod3 <- glm(rel_rec ~ region + poly(whc, 2), family = binomial(link = "logit") , data = dat ) 
mod3.1 <- glm(rel_rec ~ region + whc, family = binomial(link = "logit") , data = dat ) 
mod4 <- glm(rel_rec ~ region + poly(c_n, 2), family = binomial(link = "logit") , data = dat ) 
mod4.1 <- glm(rel_rec ~ region + c_n, family = binomial(link = "logit") , data = dat )
mod8 <- glm(rel_rec ~ region + poly(winteravgmn, 2), family = binomial(link = "logit") , data = dat ) 
mod8.1 <- glm(rel_rec ~ region + winteravgmn, family = binomial(link = "logit") , data = dat ) 
mod9 <- glm(rel_rec ~ region + poly(springtotalsnow, 2), family = binomial(link = "logit") , data = dat ) 
mod9.1 <- glm(rel_rec ~ region + springtotalsnow, family = binomial(link = "logit") , data = dat ) 

aictab(list(canopycont2 = mod1, canopycont = mod1.1, f_b2 = mod2, f_b = mod2.1, whc2 = mod3, whc = mod3.1, 
            c_n2 = mod4, c_n = mod4.1, winteravgmn2 = mod8,
            winteravgmn = mod8.1, springtotalsnow2 = mod9, springtotalsnow = mod9.1))


# Fit global model (= full model)
mod <- glm(rel_rec ~ region + canopycont + f_b + whc + c_n + winteravgmn + springtotalsnow, 
           family = binomial(link = "logit"), data = dat, na.action = "na.fail") 

# Model diagnostics
# vif(mod) #not needed with univariate models

# # CONCLUSION: no params taken out


# Variable selection with dredge
dd <- dredge(mod, evaluate = T, rank = 'AICc', extra = 'R^2', m.lim = c(0, mm))

# Look at results
subset(dd, delta <= 2)

# Save output to use in microhab_suitab.R
saveRDS(subset(dd, delta <= 2), 'outputs/candmods_cont_abilas.rds') 

avgm <- model.avg(get.models(dd, subset = delta <= 2))
saveRDS(avgm, 'outputs/avgm_cont_abilas.rds') ###*** change here


# Model diagnostics II
mod <- glm(rel_rec ~ c_n,
           family = binomial(link = "logit"), data = dat) #fit best model

plot(fitted(mod),resid(mod)) #homogeneity of variance

lag.plot(resid(mod),diag=FALSE,do.lines=FALSE) #independence of residuals to neighboring residuals

deviance(mod)/df.residual(mod) #dispersion should = 1 for variance to equal the mean

coef(lm(with(dat,rel_rec)~fitted(mod))) #link function appropriate if slope = 1

## CONCLUSION: 
# homogeneity of variance - not met
# independence of residuals - met
# variance equals the mean - not met
# link function - moderately appropriate




####################################################################################################

# # ANEOCC - MODEL AVERAGING # #
# Family: Binomial

####################################################################################################

# Data
dat <- census8 %>% 
  dplyr::filter(species == 'ANEOCC') 

# Identify how many parameters to use in model selection
mm <- round(nrow(dat)/10, 0) #rule of thumb: ~ 10 data points/parameter
mm

# Choose quadratic effects (<= 2 delta AICc)
mod1 <- glm(rel_rec ~ region + poly(canopycont, 2), family = binomial(link = "logit"), data = dat )
mod1.1 <- glm(rel_rec ~ region + canopycont, family = binomial(link = "logit") , data = dat )
mod2 <- glm(rel_rec ~ region + poly(f_b, 2), family = binomial(link = "logit") , data = dat ) 
mod2.1 <- glm(rel_rec ~ region + f_b, family = binomial(link = "logit") , data = dat ) 
mod3 <- glm(rel_rec ~ region + poly(whc, 2), family = binomial(link = "logit") , data = dat ) 
mod3.1 <- glm(rel_rec ~ region + whc, family = binomial(link = "logit") , data = dat ) 
mod4 <- glm(rel_rec ~ region + poly(c_n, 2), family = binomial(link = "logit") , data = dat ) 
mod4.1 <- glm(rel_rec ~ region + c_n, family = binomial(link = "logit") , data = dat )
mod8 <- glm(rel_rec ~ region + poly(winteravgmn, 2), family = binomial(link = "logit") , data = dat ) 
mod8.1 <- glm(rel_rec ~ region + winteravgmn, family = binomial(link = "logit") , data = dat ) 
mod9 <- glm(rel_rec ~ region + poly(springtotalsnow, 2), family = binomial(link = "logit") , data = dat ) 
mod9.1 <- glm(rel_rec ~ region + springtotalsnow, family = binomial(link = "logit") , data = dat ) 

aictab(list(canopycont2 = mod1, canopycont = mod1.1, f_b2 = mod2, f_b = mod2.1, whc2 = mod3, whc = mod3.1, 
            c_n2 = mod4, c_n = mod4.1, winteravgmn2 = mod8,
            winteravgmn = mod8.1, springtotalsnow2 = mod9, springtotalsnow = mod9.1))


# Fit global model (= full model)
mod <- glm(rel_rec ~ region + canopycont + f_b + whc + c_n + winteravgmn + springtotalsnow, 
           family = binomial(link = "logit"), data = dat, na.action = "na.fail") 

# Model diagnostics I
# vif(mod) #not needed with univariate model

# # CONCLUSION: all params used


# Variable selection with dredge
dd <- dredge(mod, evaluate = T, rank = 'AICc', extra = 'R^2', m.lim = c(0, mm))

# Look at results
subset(dd, delta <= 2)

# Save output to use in microhab_suitab.R
saveRDS(subset(dd, delta <= 2), 'outputs/candmods_cont_aneocc.rds') 

avgm <- model.avg(get.models(dd, subset = delta <= 2))
saveRDS(avgm, 'outputs/avgm_cont_aneocc.rds') ###*** change here


# Model diagnostics II
mod <- mod <- glm(rel_rec ~ f_b,
                  family = binomial(link = "logit"), data = dat) #fit best model

plot(fitted(mod),resid(mod)) #homogeneity of variance

lag.plot(resid(mod),diag=FALSE,do.lines=FALSE) #independence of residuals to neighboring residuals

deviance(mod)/df.residual(mod) #dispersion should = 1 for variance to equal the mean

coef(lm(with(dat,rel_rec)~fitted(mod))) #link function appropriate if slope = 1

## CONCLUSION: 
# homogeneity of variance - not met
# independence of residuals - not met
# variance equals the mean - not met
# link function - moderately appropriate




####################################################################################################

# # ERIPER - MODEL AVERAGING # #
# Family: Binomial

####################################################################################################

# Data
dat <- census8 %>% 
  dplyr::filter(species == 'ERIPER') 

# Identify how many parameters to use in model selection
mm <- round(nrow(dat)/10, 0) #rule of thumb: ~ 10 data points/parameter
mm

# Choose quadratic effects (<= 2 delta AICc)
mod1 <- glm(rel_rec ~ region + poly(canopycont, 2), family = binomial(link = "logit"), data = dat )
mod1.1 <- glm(rel_rec ~ region + canopycont, family = binomial(link = "logit") , data = dat )
mod2 <- glm(rel_rec ~ region + poly(f_b, 2), family = binomial(link = "logit") , data = dat ) 
mod2.1 <- glm(rel_rec ~ region + f_b, family = binomial(link = "logit") , data = dat ) 
mod3 <- glm(rel_rec ~ region + poly(whc, 2), family = binomial(link = "logit") , data = dat ) 
mod3.1 <- glm(rel_rec ~ region + whc, family = binomial(link = "logit") , data = dat ) 
mod4 <- glm(rel_rec ~ region + poly(c_n, 2), family = binomial(link = "logit") , data = dat ) 
mod4.1 <- glm(rel_rec ~ region + c_n, family = binomial(link = "logit") , data = dat )
mod8 <- glm(rel_rec ~ region + poly(winteravgmn, 2), family = binomial(link = "logit") , data = dat ) 
mod8.1 <- glm(rel_rec ~ region + winteravgmn, family = binomial(link = "logit") , data = dat ) 
mod9 <- glm(rel_rec ~ region + poly(springtotalsnow, 2), family = binomial(link = "logit") , data = dat ) 
mod9.1 <- glm(rel_rec ~ region + springtotalsnow, family = binomial(link = "logit") , data = dat ) 

aictab(list(canopycont2 = mod1, canopycont = mod1.1, f_b2 = mod2, f_b = mod2.1, whc2 = mod3, whc = mod3.1, 
            c_n2 = mod4, c_n = mod4.1, winteravgmn2 = mod8,
            winteravgmn = mod8.1, springtotalsnow2 = mod9, springtotalsnow = mod9.1))


# Fit global model (= full model)
mod <- glm(rel_rec ~ region + canopycont + poly(f_b,2) + whc + c_n + springtotalsnow, 
           family = binomial(link = "logit"), data = dat, na.action = "na.fail") 


# Model diagnostics I
vif(mod) #check for collinearity > 5

# # CONCLUSION: took out poly canopycont, winteravgm


# Variable selection with dredge
dd <- dredge(mod, evaluate = T, rank = 'AICc', extra = 'R^2', m.lim = c(0, mm))

# Look at results
subset(dd, delta <= 2)

# Save output to use in microhab_suitab.R
saveRDS(subset(dd, delta <= 2), 'outputs/candmods_cont_eriper.rds') 

avgm <- model.avg(get.models(dd, subset = delta <= 2))
saveRDS(avgm, 'outputs/avgm_cont_eriper.rds') ###*** change here


# Model diagnostics II
mod <- mod <- glm(rel_rec ~ canopycont + poly(f_b, 2),
                  family = binomial(link = "logit"), data = dat) #fit best model

plot(fitted(mod),resid(mod)) #homogeneity of variance

lag.plot(resid(mod),diag=FALSE,do.lines=FALSE) #independence of residuals to neighboring residuals

deviance(mod)/df.residual(mod) #dispersion should = 1 for variance to equal the mean

coef(lm(with(dat,rel_rec)~fitted(mod))) #link function appropriate if slope = 1

## CONCLUSION: 
# homogeneity of variance - met
# independence of residuals - met
# variance equals the mean - not met
# link function - appropriate




####################################################################################################

# # PICENG - MODEL AVERAGING # #
# Family: Binomial

####################################################################################################

# Data
dat <- census8 %>% 
  dplyr::filter(species == 'PICENG') 

# Identify how many parameters to use in model selection
mm <- round(nrow(dat)/10, 0) #rule of thumb: ~ 10 data points/parameter
mm

# Choose quadratic effects (<= 2 delta AICc)
mod1 <- glm(rel_rec ~ region + poly(canopycont, 2), family = binomial(link = "logit"), data = dat )
mod1.1 <- glm(rel_rec ~ region + canopycont, family = binomial(link = "logit") , data = dat )
mod2 <- glm(rel_rec ~ region + poly(f_b, 2), family = binomial(link = "logit") , data = dat ) 
mod2.1 <- glm(rel_rec ~ region + f_b, family = binomial(link = "logit") , data = dat ) 
mod3 <- glm(rel_rec ~ region + poly(whc, 2), family = binomial(link = "logit") , data = dat ) 
mod3.1 <- glm(rel_rec ~ region + whc, family = binomial(link = "logit") , data = dat ) 
mod4 <- glm(rel_rec ~ region + poly(c_n, 2), family = binomial(link = "logit") , data = dat ) 
mod4.1 <- glm(rel_rec ~ region + c_n, family = binomial(link = "logit") , data = dat )
mod8 <- glm(rel_rec ~ region + poly(winteravgmn, 2), family = binomial(link = "logit") , data = dat ) 
mod8.1 <- glm(rel_rec ~ region + winteravgmn, family = binomial(link = "logit") , data = dat ) 
mod9 <- glm(rel_rec ~ region + poly(springtotalsnow, 2), family = binomial(link = "logit") , data = dat ) 
mod9.1 <- glm(rel_rec ~ region + springtotalsnow, family = binomial(link = "logit") , data = dat ) 

aictab(list(canopycont2 = mod1, canopycont = mod1.1, f_b2 = mod2, f_b = mod2.1, whc2 = mod3, whc = mod3.1, 
            c_n2 = mod4, c_n = mod4.1, winteravgmn2 = mod8,
            winteravgmn = mod8.1, springtotalsnow2 = mod9, springtotalsnow = mod9.1))


# Fit global model (= full model)
mod <- glm(rel_rec ~ region + canopycont + f_b + whc + c_n + winteravgmn + springtotalsnow, 
           family = binomial(link = "logit"), data = dat, na.action = "na.fail") 

# Model diagnostics
# vif(mod) #not needed with univariate model

# # CONCLUSION: no params taken out


# Variable selection with dredge
dd <- dredge(mod, evaluate = T, rank = 'AICc', extra = 'R^2', m.lim = c(0, mm))

# Look at results
subset(dd, delta <= 2)


# Save output to use in microhab_suitab.R
saveRDS(subset(dd, delta <= 2), 'outputs/candmods_cont_piceng.rds') 

avgm <- model.avg(get.models(dd, subset = delta <= 2))
saveRDS(avgm, 'outputs/avgm_cont_piceng.rds') ###*** change here


# Model diagnostics II
mod <- glm(rel_rec ~ winteravgmn,
                  family = binomial(link = "logit"), data = dat) #fit best model

plot(fitted(mod),resid(mod)) #homogeneity of variance

lag.plot(resid(mod),diag=FALSE,do.lines=FALSE) #independence of residuals to neighboring residuals

deviance(mod)/df.residual(mod) #dispersion should = 1 for variance to equal the mean

coef(lm(with(dat,rel_rec)~fitted(mod))) #link function appropriate if slope = 1

## CONCLUSION: 
# homogeneity of variance - not met
# independence of residuals - met
# variance equals the mean - not met
# link function - appropriate




####################################################################################################

# # RUBURS - MODEL AVERAGING # #
# Family: Binomial

####################################################################################################

# Data
dat <- census8 %>% 
  dplyr::filter(species == 'RUBURS') 

# Identify how many parameters to use in model selection
mm <- round(nrow(dat)/10, 0) #rule of thumb: ~ 10 data points/parameter
mm

# Choose quadratic effects (<= 2 delta AICc)
mod1 <- glm(rel_rec ~ region + poly(canopycont, 2), family = binomial(link = "logit"), data = dat )
mod1.1 <- glm(rel_rec ~ region + canopycont, family = binomial(link = "logit") , data = dat )
mod2 <- glm(rel_rec ~ region + poly(f_b, 2), family = binomial(link = "logit") , data = dat ) 
mod2.1 <- glm(rel_rec ~ region + f_b, family = binomial(link = "logit") , data = dat ) 
mod3 <- glm(rel_rec ~ region + poly(whc, 2), family = binomial(link = "logit") , data = dat ) 
mod3.1 <- glm(rel_rec ~ region + whc, family = binomial(link = "logit") , data = dat ) 
mod4 <- glm(rel_rec ~ region + poly(c_n, 2), family = binomial(link = "logit") , data = dat ) 
mod4.1 <- glm(rel_rec ~ region + c_n, family = binomial(link = "logit") , data = dat )
mod8 <- glm(rel_rec ~ region + poly(winteravgmn, 2), family = binomial(link = "logit") , data = dat ) 
mod8.1 <- glm(rel_rec ~ region + winteravgmn, family = binomial(link = "logit") , data = dat ) 
mod9 <- glm(rel_rec ~ region + poly(springtotalsnow, 2), family = binomial(link = "logit") , data = dat ) 
mod9.1 <- glm(rel_rec ~ region + springtotalsnow, family = binomial(link = "logit") , data = dat ) 

aictab(list(canopycont2 = mod1, canopycont = mod1.1, f_b2 = mod2, f_b = mod2.1, whc2 = mod3, whc = mod3.1, 
            c_n2 = mod4, c_n = mod4.1, winteravgmn2 = mod8,
            winteravgmn = mod8.1, springtotalsnow2 = mod9, springtotalsnow = mod9.1))


# Fit global model (= full model)
mod <- glm(rel_rec ~ region + canopycont + f_b + whc + c_n + winteravgmn + springtotalsnow, 
           family = binomial(link = "logit"), data = dat, na.action = "na.fail") 


# Model diagnostics I
vif(mod) 

# # CONCLUSION: took out no params


# Variable selection with dredge
dd <- dredge(mod, evaluate = T, rank = 'AICc', extra = 'R^2', m.lim = c(0, mm))

# Look at results
subset(dd, delta <= 2)

# Save output to use in microhab_suitab.R
saveRDS(subset(dd, delta <= 2), 'outputs/candmods_cont_ruburs.rds') 

avgm <- model.avg(get.models(dd, subset = delta <= 2))
saveRDS(avgm, 'outputs/avgm_cont_ruburs.rds') ###*** change here


# Model diagnostics II
mod <- glm(rel_rec ~ whc,
           family = binomial(link = "logit"), data = dat) #fit best model

plot(fitted(mod),resid(mod)) #homogeneity of variance

lag.plot(resid(mod),diag=FALSE,do.lines=FALSE) #independence of residuals to neighboring residuals

deviance(mod)/df.residual(mod) #dispersion should = 1 for variance to equal the mean

coef(lm(with(dat,rel_rec)~fitted(mod))) #link function appropriate if slope = 1

## CONCLUSION: 
# homogeneity of variance - met
# independence of residuals - met
# variance equals the mean - not met
# link function - appropriate




####################################################################################################

# # SORSIT - MODEL AVERAGING # #
# Family: Binomial

####################################################################################################

# Data
dat <- census8 %>% 
  dplyr::filter(species == 'SORSIT') 

# Identify how many parameters to use in model selection
mm <- round(nrow(dat)/10, 0) #rule of thumb: ~ 10 data points/parameter
mm

# Choose quadratic effects (<= 2 delta AICc)
mod1 <- glm(rel_rec ~ region + poly(canopycont, 2), family = binomial(link = "logit"), data = dat )
mod1.1 <- glm(rel_rec ~ region + canopycont, family = binomial(link = "logit") , data = dat )
mod2 <- glm(rel_rec ~ region + poly(f_b, 2), family = binomial(link = "logit") , data = dat ) 
mod2.1 <- glm(rel_rec ~ region + f_b, family = binomial(link = "logit") , data = dat ) 
mod3 <- glm(rel_rec ~ region + poly(whc, 2), family = binomial(link = "logit") , data = dat ) 
mod3.1 <- glm(rel_rec ~ region + whc, family = binomial(link = "logit") , data = dat ) 
mod4 <- glm(rel_rec ~ region + poly(c_n, 2), family = binomial(link = "logit") , data = dat ) 
mod4.1 <- glm(rel_rec ~ region + c_n, family = binomial(link = "logit") , data = dat )
mod8 <- glm(rel_rec ~ region + poly(winteravgmn, 2), family = binomial(link = "logit") , data = dat ) 
mod8.1 <- glm(rel_rec ~ region + winteravgmn, family = binomial(link = "logit") , data = dat ) 
mod9 <- glm(rel_rec ~ region + poly(springtotalsnow, 2), family = binomial(link = "logit") , data = dat ) 
mod9.1 <- glm(rel_rec ~ region + springtotalsnow, family = binomial(link = "logit") , data = dat ) 

aictab(list(canopycont2 = mod1, canopycont = mod1.1, f_b2 = mod2, f_b = mod2.1, whc2 = mod3, whc = mod3.1, 
            c_n2 = mod4, c_n = mod4.1, winteravgmn2 = mod8,
            winteravgmn = mod8.1, springtotalsnow2 = mod9, springtotalsnow = mod9.1))


# Fit global model (= full model)
mod <- glm(rel_rec ~ region + poly(canopycont,2) + f_b + whc + c_n + winteravgmn + springtotalsnow, 
           family = binomial(link = "logit"), data = dat, na.action = "na.fail") 


# Model diagnostics I
vif(mod) 

# # CONCLUSION: no params removed


# Variable selection with dredge
dd <- dredge(mod, evaluate = T, rank = 'AICc', extra = 'R^2', m.lim = c(0, mm))

# Look at results
subset(dd, delta <= 2)

# Save output to use in microhab_suitab.R
saveRDS(subset(dd, delta <= 2), 'outputs/candmods_cont_sorsit.rds') 

avgm <- model.avg(get.models(dd, subset = delta <= 2))
saveRDS(avgm, 'outputs/avgm_cont_sorsit.rds') ###*** change here

# Model diagnostics II
mod <- glm(rel_rec ~ c_n,
           family = binomial(link = "logit"), data = dat) #fit best model

plot(fitted(mod),resid(mod)) #homogeneity of variance

lag.plot(resid(mod),diag=FALSE,do.lines=FALSE) #independence of residuals to neighboring residuals

deviance(mod)/df.residual(mod) #dispersion should = 1 for variance to equal the mean

coef(lm(with(dat,rel_rec)~fitted(mod))) #link function appropriate if slope = 1

## CONCLUSION: 
# homogeneity of variance - met
# independence of residuals - met
# variance equals the mean - not met
# link function - appropriate




####################################################################################################

# # TELGRA - MODEL AVERAGING # #
# Family: Binomial

####################################################################################################

# Data
dat <- census8 %>% 
  dplyr::filter(species == 'TELGRA') 

# Identify how many parameters to use in model selection
mm <- round(nrow(dat)/10, 0) #rule of thumb: ~ 10 data points/parameter
mm

# Choose quadratic effects (<= 2 delta AICc)
mod1 <- glm(rel_rec ~ region + poly(canopycont, 2), family = binomial(link = "logit"), data = dat )
mod1.1 <- glm(rel_rec ~ region + canopycont, family = binomial(link = "logit") , data = dat )
mod2 <- glm(rel_rec ~ region + poly(f_b, 2), family = binomial(link = "logit") , data = dat ) 
mod2.1 <- glm(rel_rec ~ region + f_b, family = binomial(link = "logit") , data = dat ) 
mod3 <- glm(rel_rec ~ region + poly(whc, 2), family = binomial(link = "logit") , data = dat ) 
mod3.1 <- glm(rel_rec ~ region + whc, family = binomial(link = "logit") , data = dat ) 
mod4 <- glm(rel_rec ~ region + poly(c_n, 2), family = binomial(link = "logit") , data = dat ) 
mod4.1 <- glm(rel_rec ~ region + c_n, family = binomial(link = "logit") , data = dat )
mod8 <- glm(rel_rec ~ region + poly(winteravgmn, 2), family = binomial(link = "logit") , data = dat ) 
mod8.1 <- glm(rel_rec ~ region + winteravgmn, family = binomial(link = "logit") , data = dat ) 
mod9 <- glm(rel_rec ~ region + poly(springtotalsnow, 2), family = binomial(link = "logit") , data = dat ) 
mod9.1 <- glm(rel_rec ~ region + springtotalsnow, family = binomial(link = "logit") , data = dat ) 

aictab(list(canopycont2 = mod1, canopycont = mod1.1, f_b2 = mod2, f_b = mod2.1, whc2 = mod3, whc = mod3.1, 
            c_n2 = mod4, c_n = mod4.1, winteravgmn2 = mod8,
            winteravgmn = mod8.1, springtotalsnow2 = mod9, springtotalsnow = mod9.1))


# Fit global model (= full model)
mod <- glm(rel_rec ~ region + canopycont + f_b + poly(whc, 2) + c_n + winteravgmn, 
           family = binomial(link = "logit"), data = dat, na.action = "na.fail") 


# Model diagnostics I
vif(mod) 


# # CONCLUSION: took out springtotalsnow

# Variable selection with dredge
dd <- dredge(mod, evaluate = T, rank = 'AICc', extra = 'R^2', m.lim = c(0, mm))

# Look at results
subset(dd, delta <= 2)

# Save output to use in microhab_suitab.R
saveRDS(subset(dd, delta <= 2), 'outputs/candmods_cont_telgra.rds') 

avgm <- model.avg(get.models(dd, subset = delta <= 2))
saveRDS(avgm, 'outputs/avgm_cont_telgra.rds') ###*** change here


# Model diagnostics II
mod <- glm(rel_rec ~ poly(whc,2),
           family = binomial(link = "logit"), data = dat) #fit best model

plot(fitted(mod),resid(mod)) #homogeneity of variance

lag.plot(resid(mod),diag=FALSE,do.lines=FALSE) #independence of residuals to neighboring residuals

deviance(mod)/df.residual(mod) #dispersion should = 1 for variance to equal the mean

coef(lm(with(dat,rel_rec)~fitted(mod))) #link function appropriate if slope = 1

## CONCLUSION: 
# homogeneity of variance - met
# independence of residuals - met
# variance equals the mean - not met
# link function - appropriate




####################################################################################################

# # TOLMEN - MODEL AVERAGING # #
# Family: Binomial

####################################################################################################

# Data
dat <- census8 %>% 
  dplyr::filter(species == 'TOLMEN') 

# Identify how many parameters to use in model selection
mm <- round(nrow(dat)/10, 0) #rule of thumb: ~ 10 data points/parameter
mm

# Choose quadratic effects (<= 2 delta AICc)
mod1 <- glm(rel_rec ~ region + poly(canopycont, 2), family = binomial(link = "logit"), data = dat )
mod1.1 <- glm(rel_rec ~ region + canopycont, family = binomial(link = "logit") , data = dat )
mod2 <- glm(rel_rec ~ region + poly(f_b, 2), family = binomial(link = "logit") , data = dat ) 
mod2.1 <- glm(rel_rec ~ region + f_b, family = binomial(link = "logit") , data = dat ) 
mod3 <- glm(rel_rec ~ region + poly(whc, 2), family = binomial(link = "logit") , data = dat ) 
mod3.1 <- glm(rel_rec ~ region + whc, family = binomial(link = "logit") , data = dat ) 
mod4 <- glm(rel_rec ~ region + poly(c_n, 2), family = binomial(link = "logit") , data = dat ) 
mod4.1 <- glm(rel_rec ~ region + c_n, family = binomial(link = "logit") , data = dat )
mod8 <- glm(rel_rec ~ region + poly(winteravgmn, 2), family = binomial(link = "logit") , data = dat ) 
mod8.1 <- glm(rel_rec ~ region + winteravgmn, family = binomial(link = "logit") , data = dat ) 
mod9 <- glm(rel_rec ~ region + poly(springtotalsnow, 2), family = binomial(link = "logit") , data = dat ) 
mod9.1 <- glm(rel_rec ~ region + springtotalsnow, family = binomial(link = "logit") , data = dat ) 

aictab(list(canopycont2 = mod1, canopycont = mod1.1, f_b2 = mod2, f_b = mod2.1, whc2 = mod3, whc = mod3.1, 
            c_n2 = mod4, c_n = mod4.1, winteravgmn2 = mod8,
            winteravgmn = mod8.1, springtotalsnow2 = mod9, springtotalsnow = mod9.1))


# Fit global model (= full model)
mod <- glm(rel_rec ~ region + canopycont + f_b + whc + c_n + winteravgmn, 
           family = binomial(link = "logit"), data = dat, na.action = "na.fail") 


# Model diagnostics I
vif(mod) 

# # CONCLUSION: took out springtotalsnow

# Variable selection with dredge
dd <- dredge(mod, evaluate = T, rank = 'AICc', extra = 'R^2', m.lim = c(0, mm))

# Look at results
subset(dd, delta <= 2)

# Save output to use in microhab_suitab.R
saveRDS(subset(dd, delta <= 2), 'outputs/candmods_cont_tolmen.rds') 

avgm <- model.avg(get.models(dd, subset = delta <= 2))
saveRDS(avgm, 'outputs/avgm_cont_tolmen.rds') ###*** change here


# Model diagnostics II
mod <- glm(rel_rec ~ f_b,
           family = binomial(link = "logit"), data = dat) #fit best model

plot(fitted(mod),resid(mod)) #homogeneity of variance

lag.plot(resid(mod),diag=FALSE,do.lines=FALSE) #independence of residuals to neighboring residuals

deviance(mod)/df.residual(mod) #dispersion should = 1 for variance to equal the mean

coef(lm(with(dat,rel_rec)~fitted(mod))) #link function appropriate if slope = 1

## CONCLUSION: 
# homogeneity of variance - met
# independence of residuals - met
# variance equals the mean - not met
# link function - moderately appropriate




####################################################################################################

# # VACDEL - MODEL AVERAGING # #
# Family: Binomial

####################################################################################################

# Data
dat <- census8 %>% 
  dplyr::filter(species == 'VACDEL') 

# Identify how many parameters to use in model selection
mm <- round(nrow(dat)/10, 0) #rule of thumb: ~ 10 data points/parameter
mm

# Choose quadratic effects (<= 2 delta AICc)
mod1 <- glm(rel_rec ~ region + poly(canopycont, 2), family = binomial(link = "logit"), data = dat )
mod1.1 <- glm(rel_rec ~ region + canopycont, family = binomial(link = "logit") , data = dat )
mod2 <- glm(rel_rec ~ region + poly(f_b, 2), family = binomial(link = "logit") , data = dat ) 
mod2.1 <- glm(rel_rec ~ region + f_b, family = binomial(link = "logit") , data = dat ) 
mod3 <- glm(rel_rec ~ region + poly(whc, 2), family = binomial(link = "logit") , data = dat ) 
mod3.1 <- glm(rel_rec ~ region + whc, family = binomial(link = "logit") , data = dat ) 
mod4 <- glm(rel_rec ~ region + poly(c_n, 2), family = binomial(link = "logit") , data = dat ) 
mod4.1 <- glm(rel_rec ~ region + c_n, family = binomial(link = "logit") , data = dat )
mod8 <- glm(rel_rec ~ region + poly(winteravgmn, 2), family = binomial(link = "logit") , data = dat ) 
mod8.1 <- glm(rel_rec ~ region + winteravgmn, family = binomial(link = "logit") , data = dat ) 
mod9 <- glm(rel_rec ~ region + poly(springtotalsnow, 2), family = binomial(link = "logit") , data = dat ) 
mod9.1 <- glm(rel_rec ~ region + springtotalsnow, family = binomial(link = "logit") , data = dat ) 

aictab(list(canopycont2 = mod1, canopycont = mod1.1, f_b2 = mod2, f_b = mod2.1, whc2 = mod3, whc = mod3.1, 
            c_n2 = mod4, c_n = mod4.1, winteravgmn2 = mod8,
            winteravgmn = mod8.1, springtotalsnow2 = mod9, springtotalsnow = mod9.1))


# Fit global model (= full model)
mod <- glm(rel_rec ~ region + canopycont + f_b + whc + c_n + springtotalsnow, 
           family = binomial(link = "logit"), data = dat, na.action = "na.fail") 


# Model diagnostics I
vif(mod) 

# # CONCLUSION: took out winteravgmn


# Variable selection with dredge
dd <- dredge(mod, evaluate = T, rank = 'AICc', extra = 'R^2', m.lim = c(0, mm))

# Look at results
subset(dd, delta <= 2)

# Save output to use in microhab_suitab.R
saveRDS(subset(dd, delta <= 2), 'outputs/candmods_cont_vacdel.rds') 

avgm <- model.avg(get.models(dd, subset = delta <= 2))
saveRDS(avgm, 'outputs/avgm_cont_vacdel.rds') ###*** change here


# Model diagnostics II
mod <- glm(rel_rec ~ whc,
           family = binomial(link = "logit"), data = dat) #fit best model

plot(fitted(mod),resid(mod)) #homogeneity of variance

lag.plot(resid(mod),diag=FALSE,do.lines=FALSE) #independence of residuals to neighboring residuals

deviance(mod)/df.residual(mod) #dispersion should = 1 for variance to equal the mean

coef(lm(with(dat,rel_rec)~fitted(mod))) #link function appropriate if slope = 1

## CONCLUSION: 
# homogeneity of variance - met
# independence of residuals - met
# variance equals the mean - not met
# link function - moderately appropriate




####################################################################################################

# # MAHAQU - MODEL AVERAGING # #
# Family: Binomial

####################################################################################################

# Data
dat <- census8 %>% 
  dplyr::filter(species == 'MAHAQU') 

# Identify how many parameters to use in model selection
mm <- round(nrow(dat)/10, 0) #rule of thumb: ~ 10 data points/parameter
mm

# Choose quadratic effects (<= 2 delta AICc)
mod1 <- glm(rel_rec ~ region + poly(canopycont, 2), family = binomial(link = "logit"), data = dat )
mod1.1 <- glm(rel_rec ~ region + canopycont, family = binomial(link = "logit") , data = dat )
mod2 <- glm(rel_rec ~ region + poly(f_b, 2), family = binomial(link = "logit") , data = dat ) 
mod2.1 <- glm(rel_rec ~ region + f_b, family = binomial(link = "logit") , data = dat ) 
mod3 <- glm(rel_rec ~ region + poly(whc, 2), family = binomial(link = "logit") , data = dat ) 
mod3.1 <- glm(rel_rec ~ region + whc, family = binomial(link = "logit") , data = dat ) 
mod4 <- glm(rel_rec ~ region + poly(c_n, 2), family = binomial(link = "logit") , data = dat ) 
mod4.1 <- glm(rel_rec ~ region + c_n, family = binomial(link = "logit") , data = dat )
mod8 <- glm(rel_rec ~ region + poly(winteravgmn, 2), family = binomial(link = "logit") , data = dat ) 
mod8.1 <- glm(rel_rec ~ region + winteravgmn, family = binomial(link = "logit") , data = dat ) 
mod9 <- glm(rel_rec ~ region + poly(springtotalsnow, 2), family = binomial(link = "logit") , data = dat ) 
mod9.1 <- glm(rel_rec ~ region + springtotalsnow, family = binomial(link = "logit") , data = dat ) 

aictab(list(canopycont2 = mod1, canopycont = mod1.1, f_b2 = mod2, f_b = mod2.1, whc2 = mod3, whc = mod3.1, 
            c_n2 = mod4, c_n = mod4.1, winteravgmn2 = mod8,
            winteravgmn = mod8.1, springtotalsnow2 = mod9, springtotalsnow = mod9.1))


# Fit global model (= full model)
mod <- glm(rel_rec ~ region + canopycont + f_b + whc + c_n + springtotalsnow + winteravgmn, 
           family = binomial(link = "logit"), data = dat, na.action = "na.fail") 


# Model diagnostics I
# vif(mod) #not needed with univariate model


# # CONCLUSION: no params taken out

# Variable selection with dredge
dd <- dredge(mod, evaluate = T, rank = 'AICc', extra = 'R^2', m.lim = c(0, mm))

# Look at results
subset(dd, delta <= 2)

# Save output to use in microhab_suitab.R
saveRDS(subset(dd, delta <= 2), 'outputs/candmods_cont_mahaqu.rds') 

avgm <- model.avg(get.models(dd, subset = delta <= 2))
saveRDS(avgm, 'outputs/avgm_cont_mahaqu.rds') ###*** change here


# Model diagnostics II
mod <- glm(rel_rec ~ c_n,
           family = binomial(link = "logit"), data = dat) #fit best model

plot(fitted(mod),resid(mod)) #homogeneity of variance

lag.plot(resid(mod),diag=FALSE,do.lines=FALSE) #independence of residuals to neighboring residuals

deviance(mod)/df.residual(mod) #dispersion should = 1 for variance to equal the mean

coef(lm(with(dat,rel_rec)~fitted(mod))) #link function appropriate if slope = 1

## CONCLUSION: 
# homogeneity of variance - not met
# independence of residuals - met
# variance equals the mean - not met
# link function - moderately appropriate



