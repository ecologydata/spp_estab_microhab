# Aims:
# 1. Model averaging with MuMIn for species-level models (binary recruitment data)

# Date created: 15 June 2023
# Date updated: 12 July 2023

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

# Model averaging/best model results
avgm <- read_rds('outputs/avgm_bin_mahner.rds')
candmods <- read_rds('outputs/candmods_bin_mahner.rds')

avgm <- read_rds('outputs/avgm_bin_erilan.rds') 
candmods <- read_rds('outputs/candmods_bin_erilan.rds')

avgm <- read_rds('outputs/avgm_bin_vacpar.rds')
candmods <- read_rds('outputs/candmods_bin_vacpar.rds')

avgm <- read_rds('outputs/avgm_bin_luplat.rds')
candmods <- read_rds('outputs/candmods_bin_luplat.rds')

avgm <- read_rds('outputs/avgm_bin_abilas.rds')
candmods <- read_rds('outputs/candmods_bin_abilas.rds')

avgm <- read_rds('outputs/avgm_bin_anneoc.rds')
candmods <- read_rds('outputs/candmods_bin_anneoc.rds')

avgm <- read_rds('outputs/avgm_bin_eriper.rds')
candmods <- read_rds('outputs/candmods_bin_eriper.rds')

avgm <- read_rds('outputs/avgm_bin_piceng.rds')
candmods <- read_rds('outputs/candmods_bin_piceng.rds')

avgm <- read_rds('outputs/avgm_bin_ruburs.rds')
candmods <- read_rds('outputs/candmods_bin_ruburs.rds')

avgm <- read_rds('outputs/avgm_bin_sorsit.rds')
candmods <- read_rds('outputs/candmods_bin_sorsit.rds')

avgm <- read_rds('outputs/avgm_bin_telgra.rds')
candmods <- read_rds('outputs/candmods_bin_telgra.rds')

avgm <- read_rds('outputs/avgm_bin_tolmen.rds')
candmods <- read_rds('outputs/candmods_bin_tolmen.rds')

avgm <- read_rds('outputs/avgm_bin_vacdel.rds')
candmods <- read_rds('outputs/candmods_bin_vacdel.rds')

avgm <- read_rds('outputs/avgm_bin_mahaqu.rds')
candmods <- read_rds('outputs/candmods_bin_mahaqu.rds')





####################################################################################################

# # SET UP DATA FOR MODELS (no TOMST data) # # 

####################################################################################################

# Data
load('data/census_spp.RData') #adds relative recruitment to census.tms (microhab_analysis_prep.R)

# Species with indiv in >= 8 plots
census.bin <- census.spp %>% 
  group_by(species) %>% 
  filter(species == "LUPLAT" | species == "ABILAS" | species == "ANEOCC" |
           species == "ERIPER" | species == "PICENG" |
           species == "RUBURS" | species == "SORSIT" | species == "TELGRA" |
           species == "TOLMEN" | species == "VACDEL" | species == "VACPAR" |
           species == "MAHNER" | species == "ERILAN" | species == "MAHAQU") %>% 
  mutate(recruit = if_else(new_all == 0, 0, 1)) #binary recruitment data

# Scale variables (without TOMST)
census8 <- census.bin %>%
  
  filter_at(vars(recruit, region, canopycont, f_b, whc, c_n, winteravgmn, springtotalsnow),
            all_vars(!is.na(.))) %>% #remove NAs
  
  mutate_at(vars(canopycont, f_b, whc, c_n, winteravgmn, springtotalsnow),  
            ~scale(.x)[, 1]) #turn 1 column matrix into a vector


# # Full model (with TOMST)
# mod <- glm(recruit ~ region + canopycont + f_b + whc + c_n + t1_avgmax + tmoist_avgmin + 
#             winteravgmn + springtotalsnow, family = binomial(link = "logit"), data = dat) 
# 
# 
# # Full model (without TOMST)
# mod <- lm(recruit ~ region + canopycont + f_b + whc + c_n + winteravgmn + springtotalsnow, 
#           family = binomial(link = "logit"), data = dat) 




####################################################################################################

# # MAHNER - MODEL AVERAGING (with TOMST data) # #
# Family: Binomial

####################################################################################################

# Data
dat <- census.bin %>% 
  
  filter_at(vars(recruit, region, canopycont, f_b, whc, c_n, tmoist_avgmin,
                 t1_avgmax, winteravgmn, springtotalsnow),
            all_vars(!is.na(.))) %>% #remove NAs
  
  mutate_at(vars(canopycont, f_b, whc, c_n, tmoist_avgmin,
                 t1_avgmax, winteravgmn, springtotalsnow),
             ~scale(.x)[, 1]) %>% #scale variables
  
  dplyr::filter(species == 'MAHNER') 

# Identify how many parameters to use in model selection 
mm <- round(nrow(dat)/10, 0)
mm

# Choose quadratic effects (<= 2 delta AICc)
mod1 <- glm(recruit ~ region + poly(canopycont, 2), family = binomial(link = "logit"), data = dat ) #for AICc comparison
mod1.1 <- glm(recruit ~ region + canopycont, family = binomial(link = "logit") , data = dat )
mod2 <- glm(recruit ~ region + poly(f_b, 2), family = binomial(link = "logit") , data = dat ) 
mod2.1 <- glm(recruit ~ region + f_b, family = binomial(link = "logit") , data = dat ) 
mod3 <- glm(recruit ~ region + poly(whc, 2), family = binomial(link = "logit") , data = dat ) 
mod3.1 <- glm(recruit ~ region + whc, family = binomial(link = "logit") , data = dat ) 
mod4 <- glm(recruit ~ region + poly(c_n, 2), family = binomial(link = "logit") , data = dat ) 
mod4.1 <- glm(recruit ~ region + c_n, family = binomial(link = "logit") , data = dat )
mod5 <- glm(recruit ~ region + poly(t1_avgmax, 2), family = binomial(link = "logit") , data = dat ) 
mod5.1 <- glm(recruit ~ region + t1_avgmax, family = binomial(link = "logit") , data = dat ) 
mod7 <- glm(recruit ~ region + poly(tmoist_avgmin, 2), family = binomial(link = "logit") , data = dat ) 
mod7.1 <- glm(recruit ~ region + tmoist_avgmin, family = binomial(link = "logit") , data = dat ) 
mod8 <- glm(recruit ~ region + poly(winteravgmn, 2), family = binomial(link = "logit") , data = dat ) 
mod8.1 <- glm(recruit ~ region + winteravgmn, family = binomial(link = "logit") , data = dat ) 
mod9 <- glm(recruit ~ region + poly(springtotalsnow, 2), family = binomial(link = "logit") , data = dat ) 
mod9.1 <- glm(recruit ~ region + springtotalsnow, family = binomial(link = "logit") , data = dat ) 

aictab(list(canopycont2 = mod1, canopycont = mod1.1, f_b2 = mod2, f_b = mod2.1, whc2 = mod3, whc = mod3.1, 
            c_n2 = mod4, c_n = mod4.1, t1_avgmax2 = mod5, t1_avgmax = mod5.1, 
            tmoist_avgmin2 = mod7, tmoist_avgmin = mod7.1, winteravgmn2 = mod8,
            winteravgmn = mod8.1, springtotalsnow2 = mod9, springtotalsnow = mod9.1))


# Fit global model (= full model)
mod <- glm(recruit ~ region + canopycont + f_b + whc + c_n + t1_avgmax + tmoist_avgmin + 
             winteravgmn + poly(springtotalsnow, 2), family = binomial(link = "logit") , data = dat,
           na.action = "na.fail") 


# Model diagnostics
vif(mod) 

# # CONCLUSION: kept all params

# Variable selection with dredge
dd <- dredge(mod, evaluate = T, rank = 'AICc', extra = 'R^2', m.lim = c(0, mm))

# Look at results
subset(dd, delta <= 2)

# Save output to use in microhab_suitab.R
saveRDS(subset(dd, delta <= 2), 'outputs/candmods_bin_mahner.rds') #candidate models

avgm.bin.mahner <- model.avg(get.models(dd, subset = delta <= 2)) #model averaging result
saveRDS(avgm.bin.mahner, 'outputs/avgm_bin_mahner.rds')


# Model diagnostics II
mod <- glm(recruit ~ poly(springtotalsnow, 2) + region + tmoist_avgmin + winteravgmn, 
           family = binomial(link = "logit"), data = dat) #fit best model

plot(fitted(mod),resid(mod)) #homogeneity of variance

lag.plot(resid(mod),diag=FALSE,do.lines=FALSE) #independence of residuals to neighboring residuals

deviance(mod)/df.residual(mod) #dispersion should = 1 for variance to equal the mean

coef(lm(with(dat,recruit)~fitted(mod))) #link function appropriate if slope = 1

## CONCLUSION: 
# homogeneity of variance - not met
# independence of residuals - met
# variance equals the mean - not met
# link function - not appropriate




####################################################################################################

# # ERILAN - MODEL AVERGING (with TOMST data) # #
# Family: Binomial

####################################################################################################

# Data
dat <- census.bin %>% 
  
  filter_at(vars(recruit, region, canopycont, f_b, whc, c_n, tmoist_avgmin,
                 t1_avgmax, winteravgmn, springtotalsnow),
            all_vars(!is.na(.))) %>% #remove NAs
  
  mutate_at(vars(canopycont, f_b, whc, c_n, tmoist_avgmin,
                 t1_avgmax, winteravgmn, springtotalsnow),
             ~scale(.x)[, 1]) %>% #scale variables
  
  dplyr::filter(species == 'ERILAN') 

# Identify how many parameters to use in model selection 
mm <- round(nrow(dat)/10, 0)
mm

# Choose quadratic effects (<= 2 delta AICc)
mod1 <- glm(recruit ~ region + poly(canopycont, 2), family = binomial(link = "logit"), data = dat ) #for AICc comparison
mod1.1 <- glm(recruit ~ region + canopycont, family = binomial(link = "logit") , data = dat )
mod2 <- glm(recruit ~ region + poly(f_b, 2), family = binomial(link = "logit") , data = dat ) 
mod2.1 <- glm(recruit ~ region + f_b, family = binomial(link = "logit") , data = dat ) 
mod3 <- glm(recruit ~ region + poly(whc, 2), family = binomial(link = "logit") , data = dat ) 
mod3.1 <- glm(recruit ~ region + whc, family = binomial(link = "logit") , data = dat ) 
mod4 <- glm(recruit ~ region + poly(c_n, 2), family = binomial(link = "logit") , data = dat ) 
mod4.1 <- glm(recruit ~ region + c_n, family = binomial(link = "logit") , data = dat )
mod5 <- glm(recruit ~ region + poly(t1_avgmax, 2), family = binomial(link = "logit") , data = dat ) 
mod5.1 <- glm(recruit ~ region + t1_avgmax, family = binomial(link = "logit") , data = dat ) 
mod7 <- glm(recruit ~ region + poly(tmoist_avgmin, 2), family = binomial(link = "logit") , data = dat ) 
mod7.1 <- glm(recruit ~ region + tmoist_avgmin, family = binomial(link = "logit") , data = dat ) 
mod8 <- glm(recruit ~ region + poly(winteravgmn, 2), family = binomial(link = "logit") , data = dat ) 
mod8.1 <- glm(recruit ~ region + winteravgmn, family = binomial(link = "logit") , data = dat ) 
mod9 <- glm(recruit ~ region + poly(springtotalsnow, 2), family = binomial(link = "logit") , data = dat ) 
mod9.1 <- glm(recruit ~ region + springtotalsnow, family = binomial(link = "logit") , data = dat ) 

aictab(list(canopycont2 = mod1, canopycont = mod1.1, f_b2 = mod2, f_b = mod2.1, whc2 = mod3, whc = mod3.1, 
            c_n2 = mod4, c_n = mod4.1, t1_avgmax2 = mod5, t1_avgmax = mod5.1, 
            tmoist_avgmin2 = mod7, tmoist_avgmin = mod7.1, winteravgmn2 = mod8,
            winteravgmn = mod8.1, springtotalsnow2 = mod9, springtotalsnow = mod9.1))


# Fit global model (= full model)
mod <- glm(recruit ~ region + canopycont + f_b + whc + c_n + t1_avgmax + tmoist_avgmin + 
             winteravgmn + springtotalsnow, family = binomial(link = "logit") , data = dat,
           na.action = "na.fail") 


# Model diagnostics I
vif(mod) 

# # CONCLUSION: removed poly t1_avgmax, poly springtotalsnow


# Variable selection with dredge
dd <- dredge(mod, evaluate = T, rank = 'AICc', extra = 'R^2', m.lim = c(0, mm))

# Look at results
subset(dd, delta <= 2)


# Save output to use in microhab_suitab.R
saveRDS(subset(dd, delta <= 2), 'outputs/candmods_bin_erilan.rds') #candidate models

avgm.bin.erilan <- model.avg(get.models(dd, subset = delta <= 2))
saveRDS(avgm.bin.erilan, 'outputs/avgm_bin_erilan.rds')


# Model diagnostics II
mod <- glm(recruit ~ region + whc, family = binomial(link = "logit"), data = dat) #fit best model

plot(fitted(mod),resid(mod)) #homogeneity of variance

lag.plot(resid(mod),diag=FALSE,do.lines=FALSE) #independence of residuals to neighboring residuals

deviance(mod)/df.residual(mod) #dispersion should = 1 for variance to equal the mean

coef(lm(with(dat,recruit)~fitted(mod))) #link function appropriate if slope = 1

## CONCLUSION: 
# homogeneity of variance - not met
# independence of residuals - met
# variance equals the mean - not met
# link function - not appropriate




####################################################################################################

# # VACPAR - MODEL AVERAGING # # ###*** change here
# Family: Binomial

####################################################################################################

# Data
dat <- census8 %>% 
  dplyr::filter(species == 'VACPAR') ###*** change here

# Identify how many parameters to use in model selection 
mm <- round(nrow(dat)/10, 0)
mm

# Choose quadratic effects (<= 2 delta AICc)
mod1 <- glm(recruit ~ region + poly(canopycont, 2), family = binomial(link = "logit"), data = dat ) #for AICc comparison
mod1.1 <- glm(recruit ~ region + canopycont, family = binomial(link = "logit") , data = dat )
mod2 <- glm(recruit ~ region + poly(f_b, 2), family = binomial(link = "logit") , data = dat ) 
mod2.1 <- glm(recruit ~ region + f_b, family = binomial(link = "logit") , data = dat ) 
mod3 <- glm(recruit ~ region + poly(whc, 2), family = binomial(link = "logit") , data = dat ) 
mod3.1 <- glm(recruit ~ region + whc, family = binomial(link = "logit") , data = dat ) 
mod4 <- glm(recruit ~ region + poly(c_n, 2), family = binomial(link = "logit") , data = dat ) 
mod4.1 <- glm(recruit ~ region + c_n, family = binomial(link = "logit") , data = dat )
mod8 <- glm(recruit ~ region + poly(winteravgmn, 2), family = binomial(link = "logit") , data = dat ) 
mod8.1 <- glm(recruit ~ region + winteravgmn, family = binomial(link = "logit") , data = dat ) 
mod9 <- glm(recruit ~ region + poly(springtotalsnow, 2), family = binomial(link = "logit") , data = dat ) 
mod9.1 <- glm(recruit ~ region + springtotalsnow, family = binomial(link = "logit") , data = dat ) 

aictab(list(canopycont2 = mod1, canopycont = mod1.1, f_b2 = mod2, f_b = mod2.1, whc2 = mod3, 
            whc = mod3.1, c_n2 = mod4, c_n = mod4.1, winteravgmn2 = mod8, winteravgmn = mod8.1, 
            springtotalsnow2 = mod9, springtotalsnow = mod9.1))


# Fit global model (= full model)
mod <- glm(recruit ~ region + canopycont + f_b + whc + c_n + winteravgmn + poly(springtotalsnow, 2), 
           family = binomial(link = "logit"), data = dat, na.action = "na.fail") ###*** add poly terms if identified above


# Model diagnostics I
vif(mod) 

# # CONCLUSION: removed no params ###*** change here


# Variable selection with dredge
dd <- dredge(mod, evaluate = T, rank = 'AICc', extra = 'R^2', m.lim = c(0, mm))

# Look at results
subset(dd, delta <= 2)

# Save output to use in microhab_suitab.R
saveRDS(subset(dd, delta <= 2), 'outputs/candmods_bin_vacpar.rds') #candidate models

avgm <- model.avg(get.models(dd, subset = delta <= 2))
saveRDS(avgm, 'outputs/avgm_bin_vacpar.rds') ###*** change here


# Model diagnostics II
mod <- glm(recruit ~ poly(springtotalsnow, 2) + region, ###*** change here
           family = binomial(link = "logit"), data = dat) #fit best model

plot(fitted(mod),resid(mod)) #homogeneity of variance

lag.plot(resid(mod),diag=FALSE,do.lines=FALSE) #independence of residuals to neighboring residuals

deviance(mod)/df.residual(mod) #dispersion should = 1 for variance to equal the mean

coef(lm(with(dat,recruit)~fitted(mod))) #link function appropriate if slope = 1

## CONCLUSION: 
# homogeneity of variance - not met
# independence of residuals - met
# variance equals the mean - not met
# link function - appropriate




####################################################################################################

# # LUPLAT - MODEL AVERAGING # # ###*** change here
# Family: Binomial

####################################################################################################

# Data
dat <- census8 %>% 
  dplyr::filter(species == 'LUPLAT') ###*** change here

# Identify how many parameters to use in model selection 
mm <- round(nrow(dat)/10, 0)
mm

# Choose quadratic effects (<= 2 delta AICc)
mod1 <- glm(recruit ~ region + poly(canopycont, 2), family = binomial(link = "logit"), data = dat ) #for AICc comparison
mod1.1 <- glm(recruit ~ region + canopycont, family = binomial(link = "logit") , data = dat )
mod2 <- glm(recruit ~ region + poly(f_b, 2), family = binomial(link = "logit") , data = dat ) 
mod2.1 <- glm(recruit ~ region + f_b, family = binomial(link = "logit") , data = dat ) 
mod3 <- glm(recruit ~ region + poly(whc, 2), family = binomial(link = "logit") , data = dat ) 
mod3.1 <- glm(recruit ~ region + whc, family = binomial(link = "logit") , data = dat ) 
mod4 <- glm(recruit ~ region + poly(c_n, 2), family = binomial(link = "logit") , data = dat ) 
mod4.1 <- glm(recruit ~ region + c_n, family = binomial(link = "logit") , data = dat )
mod8 <- glm(recruit ~ region + poly(winteravgmn, 2), family = binomial(link = "logit") , data = dat ) 
mod8.1 <- glm(recruit ~ region + winteravgmn, family = binomial(link = "logit") , data = dat ) 
mod9 <- glm(recruit ~ region + poly(springtotalsnow, 2), family = binomial(link = "logit") , data = dat ) 
mod9.1 <- glm(recruit ~ region + springtotalsnow, family = binomial(link = "logit") , data = dat ) 

aictab(list(canopycont2 = mod1, canopycont = mod1.1, f_b2 = mod2, f_b = mod2.1, whc2 = mod3, 
            whc = mod3.1, c_n2 = mod4, c_n = mod4.1, winteravgmn2 = mod8, winteravgmn = mod8.1, 
            springtotalsnow2 = mod9, springtotalsnow = mod9.1))


# Fit global model (= full model)
mod <- glm(recruit ~ region + canopycont + f_b + whc + c_n + winteravgmn + springtotalsnow, 
           family = binomial(link = "logit"), data = dat, na.action = "na.fail") ###*** add poly terms if identified above


# Model diagnostics I
vif(mod) 

# # CONCLUSION: removed no params ###*** change here


# Variable selection with dredge
dd <- dredge(mod, evaluate = T, rank = 'AICc', extra = 'R^2', m.lim = c(0, mm))

# Look at results
subset(dd, delta <= 2)

# Save output to use in microhab_suitab.R
saveRDS(subset(dd, delta <= 2), 'outputs/candmods_bin_luplat.rds') #candidate models

avgm <- model.avg(get.models(dd, subset = delta <= 2))
saveRDS(avgm, 'outputs/avgm_bin_luplat.rds') ###*** change here


# Model diagnostics II
mod <- glm(recruit ~ canopycont + springtotalsnow, ###*** change here
           family = binomial(link = "logit"), data = dat) #fit best model

plot(fitted(mod),resid(mod)) #homogeneity of variance

lag.plot(resid(mod),diag=FALSE,do.lines=FALSE) #independence of residuals to neighboring residuals

deviance(mod)/df.residual(mod) #dispersion should = 1 for variance to equal the mean

coef(lm(with(dat,recruit)~fitted(mod))) #link function appropriate if slope = 1

## CONCLUSION: 
# homogeneity of variance - not met
# independence of residuals - not met
# variance equals the mean - not met
# link function - moderately appropriate




####################################################################################################

# # ABILAS - MODEL AVERAGING # # ###*** change here
# Family: Binomial

####################################################################################################

# Data
dat <- census8 %>% 
  dplyr::filter(species == 'ABILAS') ###*** change here

# Identify how many parameters to use in model selection 
mm <- round(nrow(dat)/10, 0)
mm

# Choose quadratic effects (<= 2 delta AICc)
mod1 <- glm(recruit ~ region + poly(canopycont, 2), family = binomial(link = "logit"), data = dat ) #for AICc comparison
mod1.1 <- glm(recruit ~ region + canopycont, family = binomial(link = "logit") , data = dat )
mod2 <- glm(recruit ~ region + poly(f_b, 2), family = binomial(link = "logit") , data = dat ) 
mod2.1 <- glm(recruit ~ region + f_b, family = binomial(link = "logit") , data = dat ) 
mod3 <- glm(recruit ~ region + poly(whc, 2), family = binomial(link = "logit") , data = dat ) 
mod3.1 <- glm(recruit ~ region + whc, family = binomial(link = "logit") , data = dat ) 
mod4 <- glm(recruit ~ region + poly(c_n, 2), family = binomial(link = "logit") , data = dat ) 
mod4.1 <- glm(recruit ~ region + c_n, family = binomial(link = "logit") , data = dat )
mod8 <- glm(recruit ~ region + poly(winteravgmn, 2), family = binomial(link = "logit") , data = dat ) 
mod8.1 <- glm(recruit ~ region + winteravgmn, family = binomial(link = "logit") , data = dat ) 
mod9 <- glm(recruit ~ region + poly(springtotalsnow, 2), family = binomial(link = "logit") , data = dat ) 
mod9.1 <- glm(recruit ~ region + springtotalsnow, family = binomial(link = "logit") , data = dat ) 

aictab(list(canopycont2 = mod1, canopycont = mod1.1, f_b2 = mod2, f_b = mod2.1, whc2 = mod3, 
            whc = mod3.1, c_n2 = mod4, c_n = mod4.1, winteravgmn2 = mod8, winteravgmn = mod8.1, 
            springtotalsnow2 = mod9, springtotalsnow = mod9.1))


# Fit global model (= full model)
mod <- glm(recruit ~ region + canopycont + f_b + poly(whc, 2) + c_n + winteravgmn + springtotalsnow, 
           family = binomial(link = "logit"), data = dat, na.action = "na.fail") ###*** add poly terms if identified above


# Model diagnostics I
vif(mod) 


# # CONCLUSION: removed no params; ###*** change here


# Variable selection with dredge
dd <- dredge(mod, evaluate = T, rank = 'AICc', extra = 'R^2', m.lim = c(0, mm))

# Look at results
subset(dd, delta <= 2)

# Save output to use in microhab_suitab.R
saveRDS(subset(dd, delta <= 2), 'outputs/candmods_bin_abilas.rds') #candidate models

avgm <- model.avg(get.models(dd, subset = delta <= 2))
saveRDS(avgm, 'outputs/avgm_bin_abilas.rds') ###*** change here


# Model diagnostics II
mod <- glm(recruit ~ region + c_n + winteravgmn, ###*** change here
           family = binomial(link = "logit"), data = dat) #fit best model

plot(fitted(mod),resid(mod)) #homogeneity of variance

lag.plot(resid(mod),diag=FALSE,do.lines=FALSE) #independence of residuals to neighboring residuals

deviance(mod)/df.residual(mod) #dispersion should = 1 for variance to equal the mean

coef(lm(with(dat,recruit)~fitted(mod))) #link function appropriate if slope = 1

## CONCLUSION: 
# homogeneity of variance - not met
# independence of residuals - not met
# variance equals the mean - not met
# link function - appropriate




####################################################################################################

# # ANEOCC - MODEL AVERAGING # # ###*** change here
# Family: Binomial

####################################################################################################

# Data
dat <- census8 %>% 
  dplyr::filter(species == 'ANEOCC') ###*** change here

# Identify how many parameters to use in model selection 
mm <- round(nrow(dat)/10, 0)
mm

# Choose quadratic effects (<= 2 delta AICc)
mod1 <- glm(recruit ~ region + poly(canopycont, 2), family = binomial(link = "logit"), data = dat ) #for AICc comparison
mod1.1 <- glm(recruit ~ region + canopycont, family = binomial(link = "logit") , data = dat )
mod2 <- glm(recruit ~ region + poly(f_b, 2), family = binomial(link = "logit") , data = dat ) 
mod2.1 <- glm(recruit ~ region + f_b, family = binomial(link = "logit") , data = dat ) 
mod3 <- glm(recruit ~ region + poly(whc, 2), family = binomial(link = "logit") , data = dat ) 
mod3.1 <- glm(recruit ~ region + whc, family = binomial(link = "logit") , data = dat ) 
mod4 <- glm(recruit ~ region + poly(c_n, 2), family = binomial(link = "logit") , data = dat ) 
mod4.1 <- glm(recruit ~ region + c_n, family = binomial(link = "logit") , data = dat )
mod8 <- glm(recruit ~ region + poly(winteravgmn, 2), family = binomial(link = "logit") , data = dat ) 
mod8.1 <- glm(recruit ~ region + winteravgmn, family = binomial(link = "logit") , data = dat ) 
mod9 <- glm(recruit ~ region + poly(springtotalsnow, 2), family = binomial(link = "logit") , data = dat ) 
mod9.1 <- glm(recruit ~ region + springtotalsnow, family = binomial(link = "logit") , data = dat ) 

aictab(list(canopycont2 = mod1, canopycont = mod1.1, f_b2 = mod2, f_b = mod2.1, whc2 = mod3, 
            whc = mod3.1, c_n2 = mod4, c_n = mod4.1, winteravgmn2 = mod8, winteravgmn = mod8.1, 
            springtotalsnow2 = mod9, springtotalsnow = mod9.1))


# Fit global model (= full model)
mod <- glm(recruit ~ region + canopycont + f_b + whc + c_n + winteravgmn + springtotalsnow, 
           family = binomial(link = "logit"), data = dat, na.action = "na.fail") ###*** add poly terms if identified above

# Model diagnostics I
vif(mod) 

# # CONCLUSION: removed no params; ###*** change here


# Variable selection with dredge
dd <- dredge(mod, evaluate = T, rank = 'AICc', extra = 'R^2', m.lim = c(0, mm))

# Look at results
subset(dd, delta <= 2)

# Save output to use in microhab_suitab.R
saveRDS(subset(dd, delta <= 2), 'outputs/candmods_bin_anneoc.rds') #candidate models

avgm <- model.avg(get.models(dd, subset = delta <= 2))
saveRDS(avgm, 'outputs/avgm_bin_anneoc.rds') ###*** change here


# Model diagnostics II
mod <- glm(recruit ~ region + springtotalsnow, ###*** change here
           family = binomial(link = "logit"), data = dat) #fit best model

plot(fitted(mod),resid(mod)) #homogeneity of variance

lag.plot(resid(mod),diag=FALSE,do.lines=FALSE) #independence of residuals to neighboring residuals

deviance(mod)/df.residual(mod) #dispersion should = 1 for variance to equal the mean

coef(lm(with(dat,recruit)~fitted(mod))) #link function appropriate if slope = 1

## CONCLUSION: 
# homogeneity of variance - not met
# independence of residuals - not met
# variance equals the mean - not met
# link function - moderately appropriate



####################################################################################################

# # ERIPER - MODEL AVERAGING # # ###*** change here
# Family: Binomial

####################################################################################################

# Data
dat <- census8 %>% 
  dplyr::filter(species == 'ERIPER') ###*** change here

# Identify how many parameters to use in model selection
mm <- round(nrow(dat)/10, 0)
mm

# Choose quadratic effects (<= 2 delta AICc)
mod1 <- glm(recruit ~ region + poly(canopycont, 2), family = binomial(link = "logit"), data = dat ) #for AICc comparison
mod1.1 <- glm(recruit ~ region + canopycont, family = binomial(link = "logit") , data = dat )
mod2 <- glm(recruit ~ region + poly(f_b, 2), family = binomial(link = "logit") , data = dat ) 
mod2.1 <- glm(recruit ~ region + f_b, family = binomial(link = "logit") , data = dat ) 
mod3 <- glm(recruit ~ region + poly(whc, 2), family = binomial(link = "logit") , data = dat ) 
mod3.1 <- glm(recruit ~ region + whc, family = binomial(link = "logit") , data = dat ) 
mod4 <- glm(recruit ~ region + poly(c_n, 2), family = binomial(link = "logit") , data = dat ) 
mod4.1 <- glm(recruit ~ region + c_n, family = binomial(link = "logit") , data = dat )
mod8 <- glm(recruit ~ region + poly(winteravgmn, 2), family = binomial(link = "logit") , data = dat ) 
mod8.1 <- glm(recruit ~ region + winteravgmn, family = binomial(link = "logit") , data = dat ) 
mod9 <- glm(recruit ~ region + poly(springtotalsnow, 2), family = binomial(link = "logit") , data = dat ) 
mod9.1 <- glm(recruit ~ region + springtotalsnow, family = binomial(link = "logit") , data = dat ) 

aictab(list(canopycont2 = mod1, canopycont = mod1.1, f_b2 = mod2, f_b = mod2.1, whc2 = mod3, 
            whc = mod3.1, c_n2 = mod4, c_n = mod4.1, winteravgmn2 = mod8, winteravgmn = mod8.1, 
            springtotalsnow2 = mod9, springtotalsnow = mod9.1))


# Fit global model (= full model)
mod <- glm(recruit ~ region + poly(canopycont, 2) + f_b + whc + c_n + winteravgmn + springtotalsnow, 
           family = binomial(link = "logit"), data = dat, na.action = "na.fail") ###*** add poly terms if identified above

# Model diagnostics I
vif(mod) 

# # CONCLUSION: removed no params ###*** change here


# Variable selection with dredge
dd <- dredge(mod, evaluate = T, rank = 'AICc', extra = 'R^2', m.lim = c(0, mm))

# Look at results
subset(dd, delta <= 2)

# Save output to use in microhab_suitab.R
saveRDS(subset(dd, delta <= 2), 'outputs/candmods_bin_eriper.rds') #candidate models

avgm <- model.avg(get.models(dd, subset = delta <= 2))
saveRDS(avgm, 'outputs/avgm_bin_eriper.rds') ###*** change here


# Model diagnostics II
mod <- glm(recruit ~ poly(canopycont, 2) + c_n + springtotalsnow, ###*** change here
           family = binomial(link = "logit"), data = dat) #fit best model

plot(fitted(mod),resid(mod)) #homogeneity of variance

lag.plot(resid(mod),diag=FALSE,do.lines=FALSE) #independence of residuals to neighboring residuals

deviance(mod)/df.residual(mod) #dispersion should = 1 for variance to equal the mean

coef(lm(with(dat,recruit)~fitted(mod))) #link function appropriate if slope = 1

## CONCLUSION: 
# homogeneity of variance - not met
# independence of residuals - not met
# variance equals the mean - not met
# link function - appropriate




####################################################################################################

# # PICENG - MODEL AVERAGING # # ###*** change here
# Family: Binomial

####################################################################################################

# Data
dat <- census8 %>% 
  dplyr::filter(species == 'PICENG') ###*** change here

# Identify how many parameters to use in model selection
mm <- round(nrow(dat)/10, 0)
mm

# Choose quadratic effects (<= 2 delta AICc)
mod1 <- glm(recruit ~ region + poly(canopycont, 2), family = binomial(link = "logit"), data = dat ) #for AICc comparison
mod1.1 <- glm(recruit ~ region + canopycont, family = binomial(link = "logit") , data = dat )
mod2 <- glm(recruit ~ region + poly(f_b, 2), family = binomial(link = "logit") , data = dat ) 
mod2.1 <- glm(recruit ~ region + f_b, family = binomial(link = "logit") , data = dat ) 
mod3 <- glm(recruit ~ region + poly(whc, 2), family = binomial(link = "logit") , data = dat ) 
mod3.1 <- glm(recruit ~ region + whc, family = binomial(link = "logit") , data = dat ) 
mod4 <- glm(recruit ~ region + poly(c_n, 2), family = binomial(link = "logit") , data = dat ) 
mod4.1 <- glm(recruit ~ region + c_n, family = binomial(link = "logit") , data = dat )
mod8 <- glm(recruit ~ region + poly(winteravgmn, 2), family = binomial(link = "logit") , data = dat ) 
mod8.1 <- glm(recruit ~ region + winteravgmn, family = binomial(link = "logit") , data = dat ) 
mod9 <- glm(recruit ~ region + poly(springtotalsnow, 2), family = binomial(link = "logit") , data = dat ) 
mod9.1 <- glm(recruit ~ region + springtotalsnow, family = binomial(link = "logit") , data = dat ) 

aictab(list(canopycont2 = mod1, canopycont = mod1.1, f_b2 = mod2, f_b = mod2.1, whc2 = mod3, 
            whc = mod3.1, c_n2 = mod4, c_n = mod4.1, winteravgmn2 = mod8, winteravgmn = mod8.1, 
            springtotalsnow2 = mod9, springtotalsnow = mod9.1))


# Fit global model (= full model)
mod <- glm(recruit ~ region + canopycont + f_b + whc + poly(c_n, 2) + winteravgmn + 
             poly(springtotalsnow, 2),
           family = binomial(link = "logit"), data = dat, na.action = "na.fail") ###*** add poly terms if identified above


# Model diagnostics I
vif(mod) 

# # CONCLUSION: removed no params ###*** change here


# Variable selection with dredge
dd <- dredge(mod, evaluate = T, rank = 'AICc', extra = 'R^2', m.lim = c(0, mm))

# Look at results
subset(dd, delta <= 2)

# Save output to use in microhab_suitab.R
saveRDS(subset(dd, delta <= 2), 'outputs/candmods_bin_piceng.rds') #candidate models

avgm <- model.avg(get.models(dd, subset = delta <= 2))
saveRDS(avgm, 'outputs/avgm_bin_piceng.rds') ###*** change here


# Model diagnostics II
mod <- glm(recruit ~ poly(springtotalsnow, 2) + poly(c_n, 2), ###*** change here
           family = binomial(link = "logit"), data = dat) #fit best model

plot(fitted(mod),resid(mod)) #homogeneity of variance

lag.plot(resid(mod),diag=FALSE,do.lines=FALSE) #independence of residuals to neighboring residuals

deviance(mod)/df.residual(mod) #dispersion should = 1 for variance to equal the mean

coef(lm(with(dat,recruit)~fitted(mod))) #link function appropriate if slope = 1

## CONCLUSION: 
# homogeneity of variance - not met
# independence of residuals - not met
# variance equals the mean - not met
# link function - not appropriate




####################################################################################################

# # RUBURS - MODEL AVERAGING # # ###*** change here
# Family: Binomial

####################################################################################################

# Data
dat <- census8 %>% 
  dplyr::filter(species == 'RUBURS') ###*** change here

# Identify how many parameters to use in model selection 
mm <- round(nrow(dat)/10, 0)
mm

# Choose quadratic effects (<= 2 delta AICc)
mod1 <- glm(recruit ~ region + poly(canopycont, 2), family = binomial(link = "logit"), data = dat ) #for AICc comparison
mod1.1 <- glm(recruit ~ region + canopycont, family = binomial(link = "logit") , data = dat )
mod2 <- glm(recruit ~ region + poly(f_b, 2), family = binomial(link = "logit") , data = dat ) 
mod2.1 <- glm(recruit ~ region + f_b, family = binomial(link = "logit") , data = dat ) 
mod3 <- glm(recruit ~ region + poly(whc, 2), family = binomial(link = "logit") , data = dat ) 
mod3.1 <- glm(recruit ~ region + whc, family = binomial(link = "logit") , data = dat ) 
mod4 <- glm(recruit ~ region + poly(c_n, 2), family = binomial(link = "logit") , data = dat ) 
mod4.1 <- glm(recruit ~ region + c_n, family = binomial(link = "logit") , data = dat )
mod8 <- glm(recruit ~ region + poly(winteravgmn, 2), family = binomial(link = "logit") , data = dat ) 
mod8.1 <- glm(recruit ~ region + winteravgmn, family = binomial(link = "logit") , data = dat ) 
mod9 <- glm(recruit ~ region + poly(springtotalsnow, 2), family = binomial(link = "logit") , data = dat ) 
mod9.1 <- glm(recruit ~ region + springtotalsnow, family = binomial(link = "logit") , data = dat ) 

aictab(list(canopycont2 = mod1, canopycont = mod1.1, f_b2 = mod2, f_b = mod2.1, whc2 = mod3, 
            whc = mod3.1, c_n2 = mod4, c_n = mod4.1, winteravgmn2 = mod8, winteravgmn = mod8.1, 
            springtotalsnow2 = mod9, springtotalsnow = mod9.1))


# Fit global model (= full model)
mod <- glm(recruit ~ region + canopycont + f_b + whc + c_n + winteravgmn + springtotalsnow, 
           family = binomial(link = "logit"), data = dat, na.action = "na.fail") ###*** add poly terms if identified above

# Model diagnostics I
vif(mod) 

# # CONCLUSION: removed no params ###*** change here


# Variable selection with dredge
dd <- dredge(mod, evaluate = T, rank = 'AICc', extra = 'R^2', m.lim = c(0, mm))

# Look at results
subset(dd, delta <= 2)

# Save output to use in microhab_suitab.R
saveRDS(subset(dd, delta <= 2), 'outputs/candmods_bin_ruburs.rds') #candidate models

avgm <- model.avg(get.models(dd, subset = delta <= 2))
saveRDS(avgm, 'outputs/avgm_bin_ruburs.rds') ###*** change here


# Model diagnostics II
mod <- glm(recruit ~ region + winteravgmn, ###*** change here
           family = binomial(link = "logit"), data = dat) #fit best model

plot(fitted(mod),resid(mod)) #homogeneity of variance

lag.plot(resid(mod),diag=FALSE,do.lines=FALSE) #independence of residuals to neighboring residuals

deviance(mod)/df.residual(mod) #dispersion should = 1 for variance to equal the mean

coef(lm(with(dat,recruit)~fitted(mod))) #link function appropriate if slope = 1

## CONCLUSION: 
# homogeneity of variance - not met
# independence of residuals - not met
# variance equals the mean - not met
# link function - appropriate




####################################################################################################

# # SORSIT - MODEL AVERAGING # # ###*** change here
# Family: Binomial

####################################################################################################

# Data
dat <- census8 %>% 
  dplyr::filter(species == 'SORSIT') ###*** change here

# Identify how many parameters to use in model selection
mm <- round(nrow(dat)/10, 0)
mm

# Choose quadratic effects (<= 2 delta AICc)
mod1 <- glm(recruit ~ region + poly(canopycont, 2), family = binomial(link = "logit"), data = dat ) #for AICc comparison
mod1.1 <- glm(recruit ~ region + canopycont, family = binomial(link = "logit") , data = dat )
mod2 <- glm(recruit ~ region + poly(f_b, 2), family = binomial(link = "logit") , data = dat ) 
mod2.1 <- glm(recruit ~ region + f_b, family = binomial(link = "logit") , data = dat ) 
mod3 <- glm(recruit ~ region + poly(whc, 2), family = binomial(link = "logit") , data = dat ) 
mod3.1 <- glm(recruit ~ region + whc, family = binomial(link = "logit") , data = dat ) 
mod4 <- glm(recruit ~ region + poly(c_n, 2), family = binomial(link = "logit") , data = dat ) 
mod4.1 <- glm(recruit ~ region + c_n, family = binomial(link = "logit") , data = dat )
mod8 <- glm(recruit ~ region + poly(winteravgmn, 2), family = binomial(link = "logit") , data = dat ) 
mod8.1 <- glm(recruit ~ region + winteravgmn, family = binomial(link = "logit") , data = dat ) 
mod9 <- glm(recruit ~ region + poly(springtotalsnow, 2), family = binomial(link = "logit") , data = dat ) 
mod9.1 <- glm(recruit ~ region + springtotalsnow, family = binomial(link = "logit") , data = dat ) 

aictab(list(canopycont2 = mod1, canopycont = mod1.1, f_b2 = mod2, f_b = mod2.1, whc2 = mod3, 
            whc = mod3.1, c_n2 = mod4, c_n = mod4.1, winteravgmn2 = mod8, winteravgmn = mod8.1, 
            springtotalsnow2 = mod9, springtotalsnow = mod9.1))


# Fit global model (= full model)
mod <- glm(recruit ~ region + canopycont + f_b + whc + c_n + winteravgmn + springtotalsnow, 
           family = binomial(link = "logit"), data = dat, na.action = "na.fail") ###*** add poly terms if identified above


# Model diagnostics I
vif(mod) 

# # CONCLUSION: removed no params ###*** change here


# Variable selection with dredge
dd <- dredge(mod, evaluate = T, rank = 'AICc', extra = 'R^2', m.lim = c(0, mm))

# Look at results
subset(dd, delta <= 2)

# Save output to use in microhab_suitab.R
saveRDS(subset(dd, delta <= 2), 'outputs/candmods_bin_sorsit.rds') #candidate models

avgm <- model.avg(get.models(dd, subset = delta <= 2))
saveRDS(avgm, 'outputs/avgm_bin_sorsit.rds') ###*** change here


# Model diagnostics II
mod <- glm(recruit ~ springtotalsnow + region + winteravgmn, ###*** change here
           family = binomial(link = "logit"), data = dat) #fit best model

plot(fitted(mod),resid(mod)) #homogeneity of variance

lag.plot(resid(mod),diag=FALSE,do.lines=FALSE) #independence of residuals to neighboring residuals

deviance(mod)/df.residual(mod) #dispersion should = 1 for variance to equal the mean

coef(lm(with(dat,recruit)~fitted(mod))) #link function appropriate if slope = 1

## CONCLUSION: 
# homogeneity of variance - not met
# independence of residuals - not met
# variance equals the mean - not met
# link function - moderately appropriate




####################################################################################################

# # TELGRA - MODEL AVERAGING # # ###*** change here
# Family: Binomial

####################################################################################################

# Data
dat <- census8 %>% 
  dplyr::filter(species == 'TELGRA') ###*** change here

# Identify how many parameters to use in model selection 
mm <- round(nrow(dat)/10, 0)
mm

# Choose quadratic effects (<= 2 delta AICc)
mod1 <- glm(recruit ~ region + poly(canopycont, 2), family = binomial(link = "logit"), data = dat ) #for AICc comparison
mod1.1 <- glm(recruit ~ region + canopycont, family = binomial(link = "logit") , data = dat )
mod2 <- glm(recruit ~ region + poly(f_b, 2), family = binomial(link = "logit") , data = dat ) 
mod2.1 <- glm(recruit ~ region + f_b, family = binomial(link = "logit") , data = dat ) 
mod3 <- glm(recruit ~ region + poly(whc, 2), family = binomial(link = "logit") , data = dat ) 
mod3.1 <- glm(recruit ~ region + whc, family = binomial(link = "logit") , data = dat ) 
mod4 <- glm(recruit ~ region + poly(c_n, 2), family = binomial(link = "logit") , data = dat ) 
mod4.1 <- glm(recruit ~ region + c_n, family = binomial(link = "logit") , data = dat )
mod8 <- glm(recruit ~ region + poly(winteravgmn, 2), family = binomial(link = "logit") , data = dat ) 
mod8.1 <- glm(recruit ~ region + winteravgmn, family = binomial(link = "logit") , data = dat ) 
mod9 <- glm(recruit ~ region + poly(springtotalsnow, 2), family = binomial(link = "logit") , data = dat ) 
mod9.1 <- glm(recruit ~ region + springtotalsnow, family = binomial(link = "logit") , data = dat ) 

aictab(list(canopycont2 = mod1, canopycont = mod1.1, f_b2 = mod2, f_b = mod2.1, whc2 = mod3, 
            whc = mod3.1, c_n2 = mod4, c_n = mod4.1, winteravgmn2 = mod8, winteravgmn = mod8.1, 
            springtotalsnow2 = mod9, springtotalsnow = mod9.1))


# Fit global model (= full model)
mod <- glm(recruit ~ region + canopycont + f_b + whc + c_n + poly(winteravgmn, 2) + springtotalsnow, 
           family = binomial(link = "logit"), data = dat, na.action = "na.fail") ###*** add poly terms if identified above

# Model diagnostics I
vif(mod) 

# # CONCLUSION: removed no params ###*** change here


# Variable selection with dredge
dd <- dredge(mod, evaluate = T, rank = 'AICc', extra = 'R^2', m.lim = c(0, mm))

# Look at results
subset(dd, delta <= 2)

# Save output to use in microhab_suitab.R
saveRDS(subset(dd, delta <= 2), 'outputs/candmods_bin_telgra.rds') #candidate models

avgm <- model.avg(get.models(dd, subset = delta <= 2))
saveRDS(avgm, 'outputs/avgm_bin_telgra.rds') ###*** change here


# Model diagnostics II
mod <- glm(recruit ~ region + canopycont + poly(winteravgmn, 2) + c_n + f_b, ###*** change here
           family = binomial(link = "logit"), data = dat) #fit best model

plot(fitted(mod),resid(mod)) #homogeneity of variance

lag.plot(resid(mod),diag=FALSE,do.lines=FALSE) #independence of residuals to neighboring residuals

deviance(mod)/df.residual(mod) #dispersion should = 1 for variance to equal the mean

coef(lm(with(dat,recruit)~fitted(mod))) #link function appropriate if slope = 1

## CONCLUSION: 
# homogeneity of variance - not met
# independence of residuals - met
# variance equals the mean - not met
# link function - appropriate




####################################################################################################

# # TOLMEN - MODEL AVERAGING # # ###*** change here
# Family: Binomial

####################################################################################################

# Data
dat <- census8 %>% 
  dplyr::filter(species == 'TOLMEN') ###*** change here

# Identify how many parameters to use in model selection 
mm <- round(nrow(dat)/10, 0)
mm

# Choose quadratic effects (<= 2 delta AICc)
mod1 <- glm(recruit ~ region + poly(canopycont, 2), family = binomial(link = "logit"), data = dat ) #for AICc comparison
mod1.1 <- glm(recruit ~ region + canopycont, family = binomial(link = "logit") , data = dat )
mod2 <- glm(recruit ~ region + poly(f_b, 2), family = binomial(link = "logit") , data = dat ) 
mod2.1 <- glm(recruit ~ region + f_b, family = binomial(link = "logit") , data = dat ) 
mod3 <- glm(recruit ~ region + poly(whc, 2), family = binomial(link = "logit") , data = dat ) 
mod3.1 <- glm(recruit ~ region + whc, family = binomial(link = "logit") , data = dat ) 
mod4 <- glm(recruit ~ region + poly(c_n, 2), family = binomial(link = "logit") , data = dat ) 
mod4.1 <- glm(recruit ~ region + c_n, family = binomial(link = "logit") , data = dat )
mod8 <- glm(recruit ~ region + poly(winteravgmn, 2), family = binomial(link = "logit") , data = dat ) 
mod8.1 <- glm(recruit ~ region + winteravgmn, family = binomial(link = "logit") , data = dat ) 
mod9 <- glm(recruit ~ region + poly(springtotalsnow, 2), family = binomial(link = "logit") , data = dat ) 
mod9.1 <- glm(recruit ~ region + springtotalsnow, family = binomial(link = "logit") , data = dat ) 

aictab(list(canopycont2 = mod1, canopycont = mod1.1, f_b2 = mod2, f_b = mod2.1, whc2 = mod3, 
            whc = mod3.1, c_n2 = mod4, c_n = mod4.1, winteravgmn2 = mod8, winteravgmn = mod8.1, 
            springtotalsnow2 = mod9, springtotalsnow = mod9.1))


# Fit global model (= full model)
mod <- glm(recruit ~ region + canopycont + f_b + whc + c_n + winteravgmn + springtotalsnow, 
           family = binomial(link = "logit"), data = dat, na.action = "na.fail") ###*** add poly terms if identified above

# Model diagnostics I
vif(mod) 

# # CONCLUSION: removed poly(c_n) because fitted prob numerically 0 or 1 occurred & estimates seem way off ###*** change here


# Variable selection with dredge
dd <- dredge(mod, evaluate = T, rank = 'AICc', extra = 'R^2', m.lim = c(0, mm))

# Look at results
subset(dd, delta <= 2)

# Save output to use in microhab_suitab.R
saveRDS(subset(dd, delta <= 2), 'outputs/candmods_bin_tolmen.rds') #candidate models

avgm <- model.avg(get.models(dd, subset = delta <= 2))
saveRDS(avgm, 'outputs/avgm_bin_tolmen.rds') ###*** change here


# Model diagnostics II
mod <- glm(recruit ~ canopycont + c_n + springtotalsnow + whc + winteravgmn, ###*** change here
           family = binomial(link = "logit"), data = dat) #fit best model

plot(fitted(mod),resid(mod)) #homogeneity of variance

lag.plot(resid(mod),diag=FALSE,do.lines=FALSE) #independence of residuals to neighboring residuals

deviance(mod)/df.residual(mod) #dispersion should = 1 for variance to equal the mean

coef(lm(with(dat,recruit)~fitted(mod))) #link function appropriate if slope = 1

## CONCLUSION: 
# homogeneity of variance - not met
# independence of residuals - met
# variance equals the mean - not met
# link function - moderately appropriate




####################################################################################################

# # VACDEL - MODEL AVERAGING # # ###*** change here
# Family: Binomial

####################################################################################################

# Data
dat <- census8 %>% 
  dplyr::filter(species == 'VACDEL') ###*** change here

# Identify how many parameters to use in model selection 
mm <- round(nrow(dat)/10, 0)
mm

# Choose quadratic effects (<= 2 delta AICc)
mod1 <- glm(recruit ~ region + poly(canopycont, 2), family = binomial(link = "logit"), data = dat ) #for AICc comparison
mod1.1 <- glm(recruit ~ region + canopycont, family = binomial(link = "logit") , data = dat )
mod2 <- glm(recruit ~ region + poly(f_b, 2), family = binomial(link = "logit") , data = dat ) 
mod2.1 <- glm(recruit ~ region + f_b, family = binomial(link = "logit") , data = dat ) 
mod3 <- glm(recruit ~ region + poly(whc, 2), family = binomial(link = "logit") , data = dat ) 
mod3.1 <- glm(recruit ~ region + whc, family = binomial(link = "logit") , data = dat ) 
mod4 <- glm(recruit ~ region + poly(c_n, 2), family = binomial(link = "logit") , data = dat ) 
mod4.1 <- glm(recruit ~ region + c_n, family = binomial(link = "logit") , data = dat )
mod8 <- glm(recruit ~ region + poly(winteravgmn, 2), family = binomial(link = "logit") , data = dat ) 
mod8.1 <- glm(recruit ~ region + winteravgmn, family = binomial(link = "logit") , data = dat ) 
mod9 <- glm(recruit ~ region + poly(springtotalsnow, 2), family = binomial(link = "logit") , data = dat ) 
mod9.1 <- glm(recruit ~ region + springtotalsnow, family = binomial(link = "logit") , data = dat ) 

aictab(list(canopycont2 = mod1, canopycont = mod1.1, f_b2 = mod2, f_b = mod2.1, whc2 = mod3, 
            whc = mod3.1, c_n2 = mod4, c_n = mod4.1, winteravgmn2 = mod8, winteravgmn = mod8.1, 
            springtotalsnow2 = mod9, springtotalsnow = mod9.1))


# Fit global model (= full model)
mod <- glm(recruit ~ region + poly(canopycont, 2) + f_b + whc + c_n + 
             winteravgmn + springtotalsnow, 
           family = binomial(link = "logit"), data = dat, na.action = "na.fail") ###*** add poly terms if identified above


# Model diagnostics I
vif(mod) 

# # CONCLUSION: removed poly terms for c_n, f_b, winteravgmn bc fitted prob 1 or 0 occurred ###*** change here


# Variable selection with dredge
dd <- dredge(mod, evaluate = T, rank = 'AICc', extra = 'R^2', m.lim = c(0, mm))

# Look at results
subset(dd, delta <= 2)

# Save output to use in microhab_suitab.R
saveRDS(subset(dd, delta <= 2), 'outputs/candmods_bin_vacdel.rds') #candidate models

avgm <- model.avg(get.models(dd, subset = delta <= 2))
saveRDS(avgm, 'outputs/avgm_bin_vacdel.rds') ###*** change here


# Model diagnostics II
mod <- glm(recruit ~ poly(canopycont, 2) + region + winteravgmn, ###*** change here
           family = binomial(link = "logit"), data = dat) #fit best model

plot(fitted(mod),resid(mod)) #homogeneity of variance

lag.plot(resid(mod),diag=FALSE,do.lines=FALSE) #independence of residuals to neighboring residuals

deviance(mod)/df.residual(mod) #dispersion should = 1 for variance to equal the mean

coef(lm(with(dat,recruit)~fitted(mod))) #link function appropriate if slope = 1

## CONCLUSION: 
# homogeneity of variance - not met
# independence of residuals - not met
# variance equals the mean - not met
# link function - moderately appropriate




####################################################################################################

# # MAHAQU - MODEL AVERAGING # # ###*** change here
# Family: Binomial

####################################################################################################

# Data
dat <- census8 %>% 
  dplyr::filter(species == 'MAHAQU') ###*** change here

# Identify how many parameters to use in model selection 
mm <- round(nrow(dat)/10, 0)
mm

# Choose quadratic effects (<= 2 delta AICc)
mod1 <- glm(recruit ~ region + poly(canopycont, 2), family = binomial(link = "logit"), data = dat ) #for AICc comparison
mod1.1 <- glm(recruit ~ region + canopycont, family = binomial(link = "logit") , data = dat )
mod2 <- glm(recruit ~ region + poly(f_b, 2), family = binomial(link = "logit") , data = dat ) 
mod2.1 <- glm(recruit ~ region + f_b, family = binomial(link = "logit") , data = dat ) 
mod3 <- glm(recruit ~ region + poly(whc, 2), family = binomial(link = "logit") , data = dat ) 
mod3.1 <- glm(recruit ~ region + whc, family = binomial(link = "logit") , data = dat ) 
mod4 <- glm(recruit ~ region + poly(c_n, 2), family = binomial(link = "logit") , data = dat ) 
mod4.1 <- glm(recruit ~ region + c_n, family = binomial(link = "logit") , data = dat )
mod8 <- glm(recruit ~ region + poly(winteravgmn, 2), family = binomial(link = "logit") , data = dat ) 
mod8.1 <- glm(recruit ~ region + winteravgmn, family = binomial(link = "logit") , data = dat ) 
mod9 <- glm(recruit ~ region + poly(springtotalsnow, 2), family = binomial(link = "logit") , data = dat ) 
mod9.1 <- glm(recruit ~ region + springtotalsnow, family = binomial(link = "logit") , data = dat ) 

aictab(list(canopycont2 = mod1, canopycont = mod1.1, f_b2 = mod2, f_b = mod2.1, whc2 = mod3, 
            whc = mod3.1, c_n2 = mod4, c_n = mod4.1, winteravgmn2 = mod8, winteravgmn = mod8.1, 
            springtotalsnow2 = mod9, springtotalsnow = mod9.1))


# Fit global model (= full model)
mod <- glm(recruit ~ region + canopycont + f_b + whc + c_n + poly(winteravgmn, 2) + springtotalsnow, 
           family = binomial(link = "logit"), data = dat, na.action = "na.fail") ###*** add poly terms if identified above

# Model diagnostics I
vif(mod) 

# # CONCLUSION: removed no params ###*** change here


# Variable selection with dredge
dd <- dredge(mod, evaluate = T, rank = 'AICc', extra = 'R^2', m.lim = c(0, mm))

# Look at results
subset(dd, delta <= 2)

# Save output to use in microhab_suitab.R
saveRDS(subset(dd, delta <= 2), 'outputs/candmods_bin_mahaqu.rds') #candidate models

avgm <- model.avg(get.models(dd, subset = delta <= 2))
saveRDS(avgm, 'outputs/avgm_bin_mahaqu.rds') ###*** change here


# Model diagnostics II
mod <- glm(recruit ~ poly(winteravgmn, 2), ###*** change here
           family = binomial(link = "logit"), data = dat) #fit best model

plot(fitted(mod),resid(mod)) #homogeneity of variance

lag.plot(resid(mod),diag=FALSE,do.lines=FALSE) #independence of residuals to neighboring residuals

deviance(mod)/df.residual(mod) #dispersion should = 1 for variance to equal the mean

coef(lm(with(dat,recruit)~fitted(mod))) #link function appropriate if slope = 1

## CONCLUSION: 
# homogeneity of variance - not met
# independence of residuals - not met
# variance equals the mean - not met
# link function - not appropriate


