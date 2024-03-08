# Aims:
# 1. Identify biologically meaningful explanatory vars

# Date created: 15 June 2023
# Date updated: 15 June 2023

# # LIBRARIES # #
library(tidyverse)
library(dplyr)
library(lme4)
library(lmerTest)
library(MuMIn)
library(AICcmodavg)
library(car)
library(DHARMa)
library(lattice)

rm(list=ls()) 


# # INPUT FILES # #
load('data/census_spp.RData') #adds relative recruitment to census.tms (microhab_analysis_prep.R)



####################################################################################################

# # CORRELATIONS IN EXPLANATORY PARAMETERS # # 

####################################################################################################

# Check for correlation in continuous predictor variables
# Predictors with correlation > 0.7 shouldn't be used together (Tabachnick & Fidell (1996))
colnames(census.spp)
cc <- c(27, 31:34, 36:39, 57:64, 69, 71:76) #columns of predictor variables

cor(census.spp[cc], use="pairwise.complete.obs") 

# Predictor variables with corr > 0.7
# microbes ~ f_b, perc_fung, perc_bact      => use f_b
# f_b ~ perc_fung, perc_bact
# whc ~ c, h, n                             => use whc
# c ~ h, n                                  => use c_n
# t1_avgmax ~ t1_avgmin, t3_avgmax          
# t1_avgmin ~ t2_avgmax, t2_avgmin, t3_avgmax, springavgmx  => use t3_avgmin & t1_avgmax
# t2_avgmax ~ t3_avgmax
# t2_avgmin ~ t3_avgmin
# tmoist_avgmax ~ tmoist_avgmin             => use tmoist_avgmin
# winteravgmx ~ winteravgmn                 => use winteravgmn
# wintertotalsnow ~ springtotalsnow         => use springtotalsnow
# springavgmx ~ t1_avgmin, springavgmn, springotalsnow
# springavgmn ~ springtotalsnow

# Choose TOMST variables
cc <- c(57:64) #columns of TOMST variables
cor(census.spp[cc], use="pairwise.complete.obs") 

# Correlation < 0.7 for TOMST temp
# t1_avgmax ~ t2_avgmin, t3_avgmin
# t1_avgmin ~ t3_avgmin
# t2_avgmax ~ t2_avgmin, t3_avgmin
# t2_avgmin ~ t3_avgmax
# t3_avgmax ~ t3_avgmin

# Both moisture variables highly correlated => use tmoist_avgmin because drought meaningful in summer

# Choose HOBO variables
cc <- c(71:76) #columns of TOMST variables
cor(census.spp[cc], use="pairwise.complete.obs") 

# Correlation < 0.7 for HOBO vars
# winteravgmx ~ wintertotalsnow, springavgmx, springtotalsnow
# winteravgmn ~ wintertotalsnow, springavgmx, springavgmn, springtotalsnow
# wintertotalsnow ~ springavgmx, springavgmn
# springavgmx ~ 


## Conclusion: based on correlations bw predictor variables and biologically meaningful parameters => 
## region + canopycont + f_b + whc + c_n + t1_avgmax + t3_avgmin + tmoist_avgmin + winteravgmn + springtotalsnow
