# Aims:
# 1. Predict recruitment & seedling survival at each plot
# 2. Plot these predicted values for each species against elevation

# Date created: 16 May 2023
# Date updated: 12 July 2023

# # LIBRARIES # #
library(tidyverse)
library(dplyr)
library(MuMIn)
library(lattice)

rm(list=ls()) 


# # INPUT FILES # #
load('data/census_spp.RData') #adds relative recruitment to census.tms (microhab_analysis_prep.R)
# [model averaging files, see below]


# # OUTPUT FILES # #
load('outputs/pred_bin.RData') #binary recruitment predictions at plot level (microhab_suitab.R)
load('outputs/pred_cont.RData') #continuous recruitment predictions at plot level (microhab_suitab.R)
load('outputs/pred_sdlg.RData') #seedling predictions at plot level (microhab_suitab.R)




####################################################################################################

# # SUBSET TO SPECIES WITH SITES ABOVE THERMAL RANGE LIMIT # # 

####################################################################################################

# Data 
load('data/census_spp.RData') #adds relative recruitment to census.tms (microhab_analysis_prep.R)

dat <- census.spp %>% 
  filter(temprange_pos < -1)

sort(unique(dat$species)) #species with sites above thermal range limit
# ABIGRA CARSTI ERILAN LUPLAT MAHAQU MAHNER MAIDIL MAIRAC PICSIT PINPON RUBSPE
# RUBURS SAMCER TELGRA TOLMEN VACPAR

sort(unique(dat$full_species))
# "Abies grandis"         "Carex stipata"         "Eriophyllum lanatum"   "Lupinus latifolius"   
# "Mahonia aquifolium"    "Mahonia nervosa"       "Maianthemum dilatatum" "Maianthemum racemosum"
# "Picea sitchensis"      "Pinus ponderosa"       "Rubus spectabilis"     "Rubus ursinus"        
# "Sambucus caerulea"     "Tellima grandiflora"   "Tolmiea menziesii"     "Vaccinium parvifolium"


####################################################################################################

# # BINARY RECRUITMENT PREDICTIONS # # 

####################################################################################################

# Data 
load('data/census_spp.RData') #adds relative recruitment to census.tms (microhab_analysis_prep.R)

# Data for predictions with TOMST data
pred.tomst <- census.spp %>% 
  filter_at(vars(region, canopycont, f_b, whc, c_n, tmoist_avgmin,
                 t1_avgmax, winteravgmn, springtotalsnow), all_vars(!is.na(.))) %>% #remove NAs
  
  mutate_at(vars(canopycont, f_b, whc, c_n, tmoist_avgmin,
                 t1_avgmax, winteravgmn, springtotalsnow), ~scale(.x)[, 1]) #scale variables to match model


# Data for predictions without TOMST data
pred.dat <- census.spp %>% 
  filter_at(vars(region, canopycont, f_b, whc, c_n, winteravgmn, springtotalsnow), 
            all_vars(!is.na(.))) %>% 
  
  mutate_at(vars(canopycont, f_b, whc, c_n, winteravgmn, springtotalsnow), ~scale(.x)[, 1]) 

# Columns for predictions
pred.tomst$pred.recr.bin <- NA 
pred.dat$pred.recr.bin <- NA 


# Model averaging (or best model) from varsel_recruitprob_bin.R
mod <- read_rds('outputs/avgm_bin_mahner.rds') 

# Predict values for that species at each plot
pred.tomst <- pred.tomst %>% 
  mutate(pred.recr.bin = if_else(species == 'MAHNER', ###*** change here
                                 predict(mod, newdata = pred.tomst, type = 'response'), pred.recr.bin))

## PROBLEM: error when predicting with type = 'response'
## DIAGNOSIS: debugging code moved to debug_preds.R
## SOLUTION: use ~scale(.x)[, 1] to scale variables in both testing and predicting data


mod <- read_rds('outputs/avgm_bin_erilan.rds') 
pred.tomst <- pred.tomst %>% 
  mutate(pred.recr.bin = if_else(species == 'ERILAN', ###*** change here
                                 predict(mod, newdata = pred.tomst, type = 'response'), pred.recr.bin))

# Predict to non-TOMST dataframe
avgm <- read_rds('outputs/avgm_bin_vacpar.rds') 
pred.dat <- pred.dat %>% 
  mutate(pred.recr.bin = if_else(species == 'VACPAR', ###*** change here
                                 predict(avgm, newdata = pred.dat, type = 'response'), pred.recr.bin))


avgm <- read_rds('outputs/avgm_bin_luplat.rds') 
pred.dat <- pred.dat %>% 
  mutate(pred.recr.bin = if_else(species == 'LUPLAT', ###*** change here
                                 predict(avgm, newdata = pred.dat, type = 'response'), pred.recr.bin))

# avgm <- read_rds('outputs/avgm_bin_abilas.rds') 
# pred.dat <- pred.dat %>% 
#   mutate(pred.recr.bin = if_else(species == 'ABILAS', ###*** change here
#                                  predict(avgm, newdata = pred.dat), pred.recr.bin))
# 
# avgm <- read_rds('outputs/avgm_bin_anneoc.rds') 
# pred.dat <- pred.dat %>% 
#   mutate(pred.recr.bin = if_else(species == 'ANNEOC', ###*** change here
#                                  predict(avgm, newdata = pred.dat, type = 'response'), pred.recr.bin))
# 
# avgm <- read_rds('outputs/avgm_bin_eriper.rds')
# pred.dat <- pred.dat %>% 
#   mutate(pred.recr.bin = if_else(species == 'ERIPER', ###*** change here
#                                  predict(avgm, newdata = pred.dat), pred.recr.bin))
# 
# mod <- read_rds('outputs/mod_bin_piceng.rds') #best model
# summary(mod)
# nd <- filter(pred.dat, species == 'PICENG')
# nd$pred.recr.cont <- predict(mod, newdata = data.frame(springtotalsnow = nd$springtotalsnow), 
#                              type = 'response')
# 
# pred.dat <- bind_rows(pred.dat, nd) #combine with main dataframe

avgm <- read_rds('outputs/avgm_bin_ruburs.rds')
pred.dat <- pred.dat %>%
  mutate(pred.recr.bin = if_else(species == 'RUBURS', ###*** change here
                                 predict(avgm, newdata = pred.dat, type = 'response'), pred.recr.bin))

# avgm <- read_rds('outputs/avgm_bin_sorsit.rds') 
# pred.dat <- pred.dat %>% 
#   mutate(pred.recr.bin = if_else(species == 'SORSIT', ###*** change here
#                                  predict(avgm, newdata = pred.dat, type = 'response'), pred.recr.bin))

avgm <- read_rds('outputs/avgm_bin_telgra.rds') 
pred.dat <- pred.dat %>% 
  mutate(pred.recr.bin = if_else(species == 'TELGRA', ###*** change here
                                 predict(avgm, newdata = pred.dat, type = 'response'), pred.recr.bin))


avgm <- read_rds('outputs/avgm_bin_tolmen.rds') 
pred.dat <- pred.dat %>% 
  mutate(pred.recr.bin = if_else(species == 'TOLMEN', ###*** change here
                                 predict(avgm, newdata = pred.dat, type = 'response'), pred.recr.bin))

# avgm <- read_rds('outputs/avgm_bin_vacdel.rds') 
# pred.dat <- pred.dat %>% 
#   mutate(pred.recr.bin = if_else(species == 'VACDEL', ###*** change here
#                                  predict(avgm, newdata = pred.dat), pred.recr.bin))

avgm <- read_rds('outputs/avgm_bin_mahaqu.rds') 
pred.dat <- pred.dat %>% 
  mutate(pred.recr.bin = if_else(species == 'MAHAQU', ###*** change here
                                 predict(avgm, newdata = pred.dat, type = 'response'), pred.recr.bin))


# Combine dataframes
pred.bin <- bind_rows(pred.dat, pred.tomst)

pred.bin <- pred.bin %>% 
  filter(!is.na(pred.recr.bin)) #remove rows without predictions (duplicates from binding dataframes)

# Summarize results
pred.bin %>% 
  group_by(species) %>% 
  summarize(min = min(pred.recr.bin, na.rm = T), mean = mean(pred.recr.bin, na.rm = T), 
            max = max(pred.recr.bin, na.rm = T))

# Plot results
xyplot(pred.recr.bin ~ elev | species, groups = region, data = pred.bin,
       auto.key = list(draw.key = list(text = list("RP", "MB")), points = TRUE, columns = 2),
       xlab = "Elevation [m a.s.l.]",
       ylab = "Predicted recruitment probability")

# Save predictions
save(pred.bin, file = 'outputs/pred_bin.RData') #recruitment probability predictions at plot level (microhab_suitab.R)





####################################################################################################

# # CONTINUOUS RECRUITMENT PREDICTIONS # # 

####################################################################################################

# Data 
load('data/census_spp.RData') #adds relative recruitment to census.tms (microhab_analysis_prep.R)

# Data for predictions with TOMST data
pred.tomst <- census.spp %>% 
  filter_at(vars(region, canopycont, f_b, whc, c_n, tmoist_avgmin,
                 t1_avgmax, winteravgmn, springtotalsnow), all_vars(!is.na(.))) %>% #remove NAs
  
  mutate_at(vars(canopycont, f_b, whc, c_n, tmoist_avgmin,
                 t1_avgmax, winteravgmn, springtotalsnow), ~scale(.x)[, 1]) #scale variables to match model

# Data for predictions without TOMST data
pred.dat <- census.spp %>% 
  filter_at(vars(region, canopycont, f_b, whc, c_n, winteravgmn, springtotalsnow), 
            all_vars(!is.na(.))) %>% 
  
  mutate_at(vars(canopycont, f_b, whc, c_n, winteravgmn, springtotalsnow), ~scale(.x)[, 1]) 


# Set up column for predictions
pred.tomst$pred.recr.cont <- NA 
pred.dat$pred.recr.cont <- NA 


# Model averaging (or best model) from varsel_recruit_bin.R
avgm <- read_rds('outputs/avgm_cont_mahner.rds')

# Predict values for that species at each plot
pred.tomst <- pred.tomst %>% 
  mutate(pred.recr.cont = if_else(species == 'MAHNER', ###*** change here
                                  predict(avgm, newdata = pred.tomst, type = 'response'), pred.recr.cont))

avgm <- read_rds('outputs/mod_cont_erilan.rds')
pred.tomst <- pred.tomst %>% 
  mutate(pred.recr.cont = if_else(species == 'ERILAN', ###*** change here
                                  predict(avgm, newdata = pred.tomst, type = 'response'), pred.recr.cont))

# Predict to non-TOMST dataframe
avgm <- read_rds('outputs/avgm_cont_vacpar.rds')
pred.dat <- pred.dat %>% 
  mutate(pred.recr.cont = if_else(species == 'VACPAR', ###*** change here
                                  predict(avgm, newdata = pred.dat, type = 'response'), pred.recr.cont))

avgm <- read_rds('outputs/avgm_cont_luplat.rds')
pred.dat <- pred.dat %>% 
  mutate(pred.recr.cont = if_else(species == 'LUPLAT', ###*** change here
                                  predict(avgm, newdata = pred.dat, type = 'response'), pred.recr.cont))

# mod <- read_rds('outputs/mod_cont_abilas.rds') #manual entry for best model
# summary(mod)
# nd <- filter(pred.dat, species == 'ABILAS')
# nd$pred.recr.cont <- predict(mod, newdata = data.frame(c_n = nd$c_n))
# 
# pred.dat <- bind_rows(pred.dat, nd) #combine with main dataframe

# avgm <- read_rds('outputs/avgm_cont_aneocc.rds')
# pred.dat <- pred.dat %>% 
#   mutate(pred.recr.cont = if_else(species == 'ANNEOC', ###*** change here
#                                   predict(avgm, newdata = pred.dat), pred.recr.cont))

# avgm <- read_rds('outputs/avgm_cont_eriper.rds')
# pred.dat <- pred.dat %>% 
#   mutate(pred.recr.cont = if_else(species == 'ERIPER', ###*** change here
#                                   predict(avgm, newdata = pred.dat), pred.recr.cont))

# mod <- read_rds('outputs/mod_cont_piceng.rds') #best model
# summary(mod)
# nd <- filter(pred.dat, species == 'PICENG')
# nd$pred.recr.cont <- predict(mod, newdata = data.frame(winteravgmn = nd$winteravgmn))
# 
# pred.dat <- bind_rows(pred.dat, nd) #combine with main dataframe

avgm <- read_rds('outputs/avgm_cont_ruburs.rds')
pred.dat <- pred.dat %>%
  mutate(pred.recr.cont = if_else(species == 'RUBURS', ###*** change here
                                  predict(avgm, newdata = pred.dat, type = 'response'), pred.recr.cont))

# avgm <- read_rds('outputs/avgm_cont_sorsit.rds')
# pred.dat <- pred.dat %>% 
#   mutate(pred.recr.cont = if_else(species == 'SORSIT', ###*** change here
#                                   predict(avgm, newdata = pred.dat), pred.recr.cont))

mod <- read_rds('outputs/avgm_cont_telgra.rds')
pred.dat <- pred.dat %>% 
  mutate(pred.recr.cont = if_else(species == 'TELGRA', ###*** change here
                                  predict(avgm, newdata = pred.dat, type = 'response'), pred.recr.cont))

avgm <- read_rds('outputs/avgm_cont_tolmen.rds')
pred.dat <- pred.dat %>%
  mutate(pred.recr.cont = if_else(species == 'TOLMEN', ###*** change here
                                  predict(avgm, newdata = pred.dat, type = 'response'), pred.recr.cont))

# avgm <- read_rds('outputs/avgm_cont_vacdel.rds')
# pred.dat <- pred.dat %>% 
#   mutate(pred.recr.cont = if_else(species == 'VACDEL', ###*** change here
#                                   predict(avgm, newdata = pred.dat), pred.recr.cont))

avgm <- read_rds('outputs/avgm_cont_mahaqu.rds')
pred.dat <- pred.dat %>% 
  mutate(pred.recr.cont = if_else(species == 'MAHAQU', ###*** change here
                                  predict(avgm, newdata = pred.dat, type = 'response'), pred.recr.cont))

# Combine dataframes
pred.cont <- bind_rows(pred.dat, pred.tomst)

pred.cont <- pred.cont %>% 
  filter(!is.na(pred.recr.cont)) #remove rows without predictions (duplicates from binding dataframes)

# Check predictions
summary(pred.cont$pred.recr.cont) #should be btw (0,1) bc of binomial models

# Save predictions
save(pred.cont, file = 'outputs/pred_cont.RData') #continuous recruitment predictions at plot level (microhab_suitab.R)




####################################################################################################

# # SEEDLING SURVIVAL PREDICTIONS # # 

####################################################################################################

# Data 
load('data/census_sdlg.RData') #adds relative recruitment to census.tms (microhab_analysis_prep.R)

# Data for predictions with TOMST data
pred.tomst <- census.sdlg %>% 
  filter_at(vars(region, canopycont, f_b, whc, c_n, tmoist_avgmin,
                 t3_avgmin, winteravgmn, springtotalsnow), all_vars(!is.na(.))) %>% #remove NAs
  
  mutate_at(vars(canopycont, f_b, whc, c_n, tmoist_avgmin,
                 t3_avgmin, winteravgmn, springtotalsnow), ~scale(.x)[, 1]) #scale variables to match model

# Data for predictions without TOMST data
pred.dat <- census.sdlg %>% 
  filter_at(vars(region, canopycont, f_b, whc, c_n, winteravgmn, springtotalsnow), 
            all_vars(!is.na(.))) %>% 
  
  mutate_at(vars(canopycont, f_b, whc, c_n, winteravgmn, springtotalsnow), ~scale(.x)[, 1]) 

# Columns for predictions
pred.tomst$pred.1yr.surv <- NA 
pred.dat$pred.1yr.surv <- NA 

pred.tomst$pred.2yr.surv <- NA 
pred.dat$pred.2yr.surv <- NA 


## Predict 1-year seedling survival 

# Model averaging (or best model) from MuMIn_varsel.R
avgm <- read_rds('outputs/avgm_1yr_ruburs.rds') 

# Predict values for that species at each plot
pred.tomst <- pred.tomst %>% 
  mutate(pred.1yr.surv = if_else(species == 'RUBURS', ###*** change here
                                 predict(avgm, newdata = pred.tomst, type = 'response'), pred.1yr.surv))

avgm <- read_rds('outputs/avgm_1yr_telgra.rds') 
pred.dat <- pred.dat %>% 
  mutate(pred.1yr.surv = if_else(species == 'TELGRA', ###*** change here
                                 predict(avgm, newdata = pred.dat, type = 'response'), pred.1yr.surv))

avgm <- read_rds('outputs/avgm_1yr_tolmen.rds')
pred.dat <- pred.dat %>%
  mutate(pred.1yr.surv = if_else(species == 'TOLMEN', ###*** change here
                                 predict(avgm, newdata = pred.dat, type = 'response'), pred.1yr.surv))

avgm <- read_rds('outputs/avgm_1yr_luplat.rds') 
pred.dat <- pred.dat %>% 
  mutate(pred.1yr.surv = if_else(species == 'LUPLAT', ###*** change here
                                 predict(avgm, newdata = pred.dat, type = 'response'), pred.1yr.surv))

avgm <- read_rds('outputs/avgm_1yr_vacpar.rds') 
pred.dat <- pred.dat %>% 
  mutate(pred.1yr.surv = if_else(species == 'VACPAR', ###*** change here
                                 predict(avgm, newdata = pred.dat, type = 'response'), pred.1yr.surv))


## Predict 2-year seedling survival 

avgm <- read_rds('outputs/mod_2yr_ruburs.rds') 
pred.tomst <- pred.tomst %>% 
  mutate(pred.2yr.surv = if_else(species == 'RUBURS', ###*** change here
                                 predict(avgm, newdata = pred.tomst, type = 'response'), pred.2yr.surv))

avgm <- read_rds('outputs/avgm_2yr_telgra.rds') 
pred.dat <- pred.dat %>% 
  mutate(pred.2yr.surv = if_else(species == 'TELGRA', ###*** change here
                                 predict(avgm, newdata = pred.dat, type = 'response'), pred.2yr.surv))

mod <- read_rds('outputs/mod_2yr_tolmen.rds') #best model
pred.dat <- pred.dat %>%
  mutate(pred.2yr.surv = if_else(species == 'TOLMEN', ###*** change here
                                 predict(mod, newdata = pred.dat, type = 'response'), pred.2yr.surv))

avgm <- read_rds('outputs/avgm_2yr_luplat.rds') 
pred.dat <- pred.dat %>% 
  mutate(pred.2yr.surv = if_else(species == 'LUPLAT', ###*** change here
                                 predict(avgm, newdata = pred.dat, type = 'response'), pred.2yr.surv))

# avgm <- read_rds('outputs/avgm_1yr_sorsit.rds') 
# pred.dat <- pred.dat %>% 
#   mutate(pred.1yr.surv = if_else(species == 'SORSIT', ###*** change here
#                                  predict(avgm, newdata = pred.dat), pred.1yr.surv), type = 'response')
# 
# avgm <- read_rds('outputs/avgm_2yr_sorsit.rds') 
# pred.dat <- pred.dat %>% 
#   mutate(pred.2yr.surv = if_else(species == 'SORSIT', ###*** change here
#                                  predict(avgm, newdata = pred.dat), pred.2yr.surv), type = 'response')
#
# avgm <- read_rds('outputs/avgm_1yr_eriper.rds') 
# pred.dat <- pred.dat %>% 
#   mutate(pred.1yr.surv = if_else(species == 'ERIPER', ###*** change here
#                                  predict(avgm, newdata = pred.dat), pred.1yr.surv), type = 'response')
# 
# mod <- read_rds('outputs/mod_2yr_eriper.rds') #best model
# pred.dat <- pred.dat %>% 
#   mutate(pred.2yr.surv = if_else(species == 'ERIPER', ###*** change here
#                                  predict(avgm, newdata = pred.dat), pred.2yr.surv), type = 'response')
#
# avgm <- read_rds('outputs/avgm_1yr_vacdel.rds') 
# pred.dat <- pred.dat %>% 
#   mutate(pred.1yr.surv = if_else(species == 'VACDEL', ###*** change here
#                                  predict(avgm, newdata = pred.dat), pred.1yr.surv), type = 'response')
# 

# Combine dataframes
pred.sdlg <- bind_rows(pred.dat, pred.tomst)

pred.sdlg <- pred.sdlg %>% 
  filter(!is.na(pred.1yr.surv)) #remove rows without predictions (duplicates from binding dataframes)

# Check predictions
summary(pred.sdlg$pred.1yr.surv) #should be binary bc of binomial models
summary(pred.sdlg$pred.2yr.surv) 

# Save predictions
save(pred.sdlg, file = 'outputs/pred_sdlg.RData') #continuous recruitment predictions at plot level (microhab_suitab.R)


# Plot results
pred.sdlg.1yr <- pred.sdlg %>% 
  filter(!is.na(pred.1yr.surv)) #remove rows without predictions (duplicates from binding dataframes)

xyplot(pred.1yr.surv ~ elev | species, groups = region, data = pred.sdlg.1yr,
       auto.key = list(draw.key = list(text = list("RP", "MB")), points = TRUE, columns = 2),
       xlab = "Elevation [m a.s.l.]",
       ylab = "Predicted proportion of seedlings surviving 1 year")

pred.sdlg.2yr <- pred.sdlg %>% 
  filter(!is.na(pred.2yr.surv)) #remove rows without predictions (duplicates from binding dataframes)

xyplot(pred.2yr.surv ~ elev | species, groups = region, data = pred.sdlg.2yr,
       auto.key = list(draw.key = list(text = list("RP", "MB")), points = TRUE, columns = 2),
       xlab = "Elevation [m a.s.l.]",
       ylab = "Predicted proportion of seedlings surviving 2 years")





