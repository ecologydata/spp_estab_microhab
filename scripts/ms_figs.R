# Aims: MANUSCRIPT FIGURES

# Author: Nathalie Chardon, Kavya Pradhan
# Date created: 11 May 2023
# Date updated: 27 Feb 2024


# # LIBRARIES # #
library(tidyverse)
library(dplyr)
library(data.table)
library(readxl)
library(ggplot2)
library(grid)
library(gridExtra)
library(gridtext)
library(ggpubr)
library(lmerTest)
library(ggeffects)
library(ggfortify)
library(stats)


rm(list=ls()) 

# # WORKING DIRECTORIES # #
# indir <- "/Users/kavyapradhan/Documents/GitHub/CL_seedaddition/"
# setwd('/Users/kavyapradhan/Documents/GitHub/CL_seedaddition/')


# # INPUT FILES # #
load('outputs/mod_avg_June2023/ALL_mod_avg.RData') #model averaging results for all models (ms_tables.R)
load('data/tidy/census_spp.RData') #adds relative recruitment to census.tms (microhab_analysis_prep.R)
load('data/tidy/thermal2elevRL.RData') #thermal to elev range limit per species (therm2elevRL.R)
load('outputs/ms_results_Jul2023/pred_bin.RData') #binary recruitment predictions at plot level (microhab_suitab.R)
load('outputs/ms_results_Jul2023/pred_cont.RData') #continuous recruitment predictions at plot level (microhab_suitab.R)
load('outputs/ms_results_Jul2023/pred_sdlg.RData') #seedling predictions at plot level (microhab_suitab.R)


# # OUTPUT FILES # #
# [figure files, see below]



####################################################################################################

# # FIG. 1 [CONCEPTUAL] # # 

####################################################################################################

# [created in PP on NC's computer: Cascades/manuscript/figs/conceptual.pptx]

####################################################################################################

# # FIG. 2 [MICROHAB-ELEV] # # 

####################################################################################################

# Data
load('data/tidy/census_spp.RData') #adds relative recruitment to census.tms (microhab_analysis_prep.R)

dat <- census.spp %>% 
  distinct(census.spp$rep, .keep_all = TRUE) %>% #keep only rep-level data
  select(c(region, site, site1, rep, elev, MAT, MAP, canopycont, c_n, f_b, springtotalsnow, tmoist_avgmin, 
           t1_avgmax, t3_avgmin, whc, winteravgmn)) #keep only microhab vars

# No. of loggers
foo <- dat %>% 
  distinct(dat$site, .keep_all = T) #remove duplicate values in 1 block

nrow(foo[!is.na(foo$tmoist_avgmin),])
nrow(foo[!is.na(foo$t1_avgmax),])
nrow(foo[!is.na(foo$t3_avgmin),])
nrow(foo[!is.na(foo$winteravgmn),])
nrow(foo[!is.na(foo$springtotalsnow),])


# Theme
pp <- 2 #point size
tt <- 12 #text size
cc <- c('#008000', 'black') #points color (https://html-color.codes/green)

mytheme <- theme(axis.title.x = element_text(size = tt), axis.text.x = element_text(colour = 'black', size = tt),
                 axis.title.y = element_text(size = tt), axis.text.y = element_text(colour = 'black', size = tt),
                 plot.title = element_text(size = tt),
                 plot.title.position = "plot", #move title to left-hand plot space
                 panel.background = element_blank(), axis.line = element_line(colour = "black"),
                 legend.position = "none"
) 


# # Climate-unrelated microhabitat variables

# CANOPYCONT

# Calculate variance partition coefficient = variance at a given level of the model, 
# divided by the total variance (the sum of the variance parameters)
mod <- lmer(canopycont ~ 1 + (1|site1), data = dat)
vc <- VarCorr(mod) %>%
  as_tibble() %>%
  mutate(icc = vcov/sum(vcov)) %>% #variance/total variance
  select(grp, icc)

# Test microhab ~ elev relationship
mod <- lmer(canopycont ~ poly(elev, 2) * region + (1|site1), data = dat) #effect of elev
summary(mod) #only add fitted line if elev is significant, add two lines if region*elev is sign

# Plot results
title <- bquote('(a)' ~
                  # sigma['within']^2 == .(round(vc[1,2], 2)$icc)*',' ~ #remove for (1|site) only models
                  sigma['among']^2 == .(round(vc[2,2], 2)$icc))

canopycont <- dat %>%
  ggplot(aes(x = elev, y = canopycont, color = region)) + 
  geom_point(size = pp) + scale_color_manual(values = cc) +
  labs(y = "Canopy Openness\n [% cover]", x = NULL, title = title, color = 'Transect') +
  geom_smooth(method = "lm", se = FALSE) +
  mytheme


# F_B
# Variance partitioning
mod <- lmer(f_b ~ 1 + (1|site1/site), data = dat) 
vc <- VarCorr(mod) %>%
  as_tibble() %>%
  mutate(icc = vcov/sum(vcov)) %>% 
  select(grp, icc)

# Test microhab ~ elev relationship
mod <- lmer(f_b ~ poly(elev, 2) * region + (1|site1), data = dat) 
summary(mod)

# Plot results
title <- bquote('(b)' ~
                  sigma['within']^2 == .(round(vc[1,2], 2)$icc)*',' ~ #remove for (1|site) only models
                  sigma['among']^2 == .(round(vc[2,2], 2)$icc))

f_b <- dat %>% 
  ggplot(aes(x = elev, y = f_b, color = region)) + 
  geom_point(size = pp) + scale_color_manual(values = cc) +
  labs(y = "Soil Fungus:Bacteria\n", x = NULL, title = title) +
  mytheme + theme(legend.position = "none", plot.title = element_text(face = 'bold'))


# C_N
# Variance partitioning
mod <- lmer(c_n ~ 1 + (1|site1/site), data = dat) 
vc <- VarCorr(mod) %>%
  as_tibble() %>%
  mutate(icc = vcov/sum(vcov)) %>% 
  select(grp, icc)

# Test microhab ~ elev relationship
mod <- lmer(c_n ~ poly(elev, 2) * region + (1|site1), data = dat) 
summary(mod)

# Plot results
title <- bquote('(c)' ~
                  sigma['within']^2 == .(round(vc[1,2], 2)$icc)*',' ~ #remove for (1|site) only models
                  sigma['among']^2 == .(round(vc[2,2], 2)$icc))

c_n <- dat %>%
  ggplot(aes(x = elev, y = c_n, color = region)) + 
  geom_point(size = pp) + scale_color_manual(values = cc) +
  labs(y = "Soil Carbon:Nitrogen\n", x = NULL, title = title) +
  mytheme + theme(legend.position = "none")


# WHC
# Variance partitioning
mod <- lmer(whc ~ 1 + (1|site1/site), data = dat) 
vc <- VarCorr(mod) %>%
  as_tibble() %>%
  mutate(icc = vcov/sum(vcov)) %>% 
  select(grp, icc)

# Test microhab ~ elev relationship
mod <- lmer(whc ~ poly(elev, 2) * region + (1|site1), data = dat) 
summary(mod)

# Plot results
title <- bquote('(d)' ~
                  sigma['within']^2 == .(round(vc[1,2], 2)$icc)*',' ~ #remove for (1|site) only models
                  sigma['among']^2 == .(round(vc[2,2], 2)$icc))

whc <- dat %>%
  ggplot(aes(x = elev, y = whc, color = region)) + 
  geom_point(size = pp) + scale_color_manual(values = cc) +
  labs(y = "Water Holding Capacity\n", x = NULL, title = title) +
  mytheme + theme(legend.position = "none")


# Climate-related microhabitat variables

# TMOIST_AVGMIN
# Variance partitioning
mod <- lmer(tmoist_avgmin ~ 1 + (1|site1), data = dat) 
vc <- VarCorr(mod) %>%
  as_tibble() %>%
  mutate(icc = vcov/sum(vcov)) %>% 
  select(grp, icc)

# Test microhab ~ elev relationship
mod <- lmer(tmoist_avgmin ~ poly(elev, 2) * region + (1|site1), data = dat) 
summary(mod)

# Plot results
title <- bquote('(e)' ~
                  # sigma['within']^2 == .(round(vc[1,2], 2)$icc)*',' ~ #remove for (1|site) only models
                  sigma['among']^2 == .(round(vc[1,2], 2)$icc))

tmoist_avgmin <- dat %>%
  ggplot(aes(x = elev, y = tmoist_avgmin, color = region)) +
  geom_point(size = pp) + scale_color_manual(values = cc) +
  labs(y = "Summer Min. Soil\n Moisture", 
       x = NULL, title = title) +
  mytheme + theme(legend.position = "none")


# T1_AVGMAX
# Variance partitioning
mod <- lmer(t1_avgmax ~ 1 + (1|site1), data = dat) 
vc <- VarCorr(mod) %>%
  as_tibble() %>%
  mutate(icc = vcov/sum(vcov)) %>% 
  select(grp, icc)

# Test microhab ~ elev relationship
mod <- lmer(t1_avgmax ~ poly(elev, 2) * region + (1|site1), data = dat) 
summary(mod)

# Plot results
title <- bquote('(f)' ~
                  sigma['within']^2 == .(round(vc[1,2], 2)$icc)*',' ~ #remove for (1|site) only models
                  sigma['among']^2 == .(round(vc[2,2], 2)$icc))

t1_avgmax <- dat %>%
  ggplot(aes(x = elev, y = t1_avgmax, color = region)) +
  geom_point(size = pp) + scale_color_manual(values = cc) +
  labs(y = "Summer Max.\n Soil Temp [ºC]", x = NULL, title = title) +
  mytheme + theme(legend.position = "none", plot.title = element_text(face = 'bold'))


# T3_AVGMIN
# Variance partitioning
mod <- lmer(t3_avgmin ~ 1 + (1|site1), data = dat) 
vc <- VarCorr(mod) %>%
  as_tibble() %>%
  mutate(icc = vcov/sum(vcov)) %>% 
  select(grp, icc)

# Test microhab ~ elev relationship
mod <- lmer(t3_avgmin ~ poly(elev, 2) * region + (1|site1), data = dat) 
summary(mod)

# Plot results
title <- bquote('(g)' ~
                  # sigma['within']^2 == .(round(vc[1,2], 2)$icc)*',' ~ #remove for (1|site) only models
                  sigma['among']^2 == .(round(vc[1,2], 2)$icc))

t3_avgmin <- dat %>%
  ggplot(aes(x = elev, y = t3_avgmin, color = region)) +
  geom_point(size = pp) + scale_color_manual(values = cc) +
  labs(y = "Summer Min.\n Plant-height Temp [ºC]", x = NULL, title = title) +
  mytheme + theme(legend.position = "none")


# WINTERAVGMN
# Variance partitioning
mod <- lmer(winteravgmn ~ 1 + (1|site1), data = dat) 
vc <- VarCorr(mod) %>%
  as_tibble() %>%
  mutate(icc = vcov/sum(vcov)) %>% 
  select(grp, icc)

# Test microhab ~ elev relationship
mod <- lmer(winteravgmn ~ poly(elev, 2) * region + (1|site1), data = dat) 
summary(mod)

# Plot results
title <- bquote('(h)' ~
                  # sigma['within']^2 == .(round(vc[1,2], 2)$icc)*',' ~ #remove for (1|site) only models
                  sigma['among']^2 == .(round(vc[1,2], 2)$icc))

winteravgmn <- dat %>%
  ggplot(aes(x = elev, y = winteravgmn, color = region)) + 
  geom_point(size = pp) + scale_color_manual(values = cc) +
  labs(y = "Winter Min.\n Soil Temp [ºC]", x = NULL, title = title) + 
  geom_smooth(method = "lm", se = FALSE) +
  mytheme + theme(legend.position = "none")


# SPRINGTOTALSNOW
# Variance partitioning
mod <- lmer(springtotalsnow ~ 1 + (1|site1), data = dat) 
vc <- VarCorr(mod) %>%
  as_tibble() %>%
  mutate(icc = vcov/sum(vcov)) %>% 
  select(grp, icc)

# Test microhab ~ elev relationship
mod <- lmer(springtotalsnow ~ poly(elev, 2) * region + (1|site1), data = dat) 
summary(mod)

# Plot results
title <- bquote('(i)' ~
                  # sigma['within']^2 == .(round(vc[1,2], 2)$icc)*',' ~ #remove for (1|site) only models
                  sigma['among']^2 == .(round(vc[1,2], 2)$icc))

springtotalsnow <- dat %>%
  ggplot(aes(x = elev, y = springtotalsnow, color = region)) + 
  geom_point(size = pp) + scale_color_manual(values = cc, labels = c("West", "East")) +
  labs(y = "Spring Snow\n [days]", x = NULL, title = title) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, data = subset(dat, region == "MB")) + 
  mytheme + theme(plot.title = element_text(face = 'italic'), legend.position = 'bottom') + guides(color=guide_legend(title = 'Transect', title.position = 'top', nrow=1, byrow=TRUE)) 

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

plotLegend <- g_legend(springtotalsnow)

# Combine plots and save
fig <- grid.arrange(arrangeGrob(canopycont, f_b, c_n, whc, tmoist_avgmin, t1_avgmax, t3_avgmin, 
                    winteravgmn, springtotalsnow + theme(legend.position = "none"), nrow = 3, heights = c(3,3,3)), plotLegend, heights = c(10,1),
                     bottom = textGrob('Elevation [m]', gp = gpar(fontsize = tt))
)

#ggsave('outputs/ms_results_Jul2023/microhab_quadelev.jpeg', fig, width = 8, height = 8, units = c('in'))
ggsave('outputs/ms_revisions_Dec2023/microhab_quadelev.pdf', fig, width = 8, height = 8, units = c('in'))



####################################################################################################

# # FIG. 3 [PARAM-SUMMARIES] # # 

####################################################################################################

# Data
load('outputs/ms_results_Jul2023/ALL_mod_avg_SIMPLE.RData') #model averaging results for all models (ms_tables.R)

# Number of times each parameter used in model in long DF format 
#making dd.tab.simple df long so that parameter values can be grouped by +/- in stacked bar graph
dd.simple.long <- gather(dd.simple, parameter, effect, 'Year':'Spring Days with Snow^2', factor_key=TRUE) 

is.na(dd.simple.long$effect) <- dd.simple.long$effect == 0 #make 0 values into NAs

dd.simple.long <- drop_na(dd.simple.long) #drop NA values from df


# Change parameter names to be shorter in figures
params <- unique(dd.simple.long$parameter)
dd.simple.long <- dd.simple.long %>%
  mutate(parameter = if_else(parameter == as.character(params[3]), 'Canopy', parameter)) %>%
  mutate(parameter = if_else(parameter == as.character(params[4]), 'Canopy^2', parameter)) %>%
  mutate(parameter = if_else(parameter == as.character(params[5]), 'F:B', parameter)) %>%
  mutate(parameter = if_else(parameter == as.character(params[6]), 'F:B^2', parameter)) %>%
  mutate(parameter = if_else(parameter == as.character(params[7]), 'C:N', parameter)) %>%
  mutate(parameter = if_else(parameter == as.character(params[8]), 'C:N^2', parameter)) %>%
  mutate(parameter = if_else(parameter == as.character(params[9]), 'WHC', parameter)) %>%
  mutate(parameter = if_else(parameter == as.character(params[10]), 'WHC^2', parameter)) %>%
  mutate(parameter = if_else(parameter == as.character(params[11]), 'S Temp', parameter)) %>%
  mutate(parameter = if_else(parameter == as.character(params[12]), 'S Soil Moist', parameter)) %>%
  mutate(parameter = if_else(parameter == as.character(params[13]), 'W Soil Temp', parameter)) %>%
  mutate(parameter = if_else(parameter == as.character(params[14]), 'W Soil Temp^2', parameter)) %>%
  mutate(parameter = if_else(parameter == as.character(params[15]), 'Snow', parameter)) %>%
  mutate(parameter = if_else(parameter == as.character(params[16]), 'Snow^2', parameter))

# Vector of models only used for 2 spp in recruitprob/recruit counts and 1 spp in seedling surv prob
pp <- c('S Temp', 'S Soil Moist')

# Number of params used in models for each life stage (n needed to configure plots below)
dd.simple %>% count(Model)

# Theme
tt <- 16
tt.x <- 12
mytheme <- theme(axis.text.x = element_text(colour = 'black', size = tt.x, angle = 60, margin = margin(t = 10), hjust = 1),
                 axis.title.y = element_text(size = tt), axis.text.y = element_text(colour = 'black', size = tt),
                 plot.title = element_text(size = tt),
                 plot.tag = element_text(size = tt), plot.tag.position = c(0.1, 0.97),
                 panel.background = element_blank(), axis.line = element_line(colour = "black"),
                 legend.text = element_text(size = tt), legend.title = element_blank(),
                 legend.key = element_rect(fill = "transparent"),
                 legend.position = c(0.001, 0.999), legend.justification = c(0, 1), 
                 legend.margin = margin(0, 6, 3, 0), legend.box.background = element_rect(colour = 'black')
) 

cc2 <- c('#AA4499', '#008000') #color-blind friendly palette to use for +/- effect 


# # Summary plots

# Recruitment probability

n <- dd.simple.long %>% 
  filter(Model == 'recruitment probability') %>% 
  distinct(Species) %>% 
  nrow() #number of models for title

recprob <- dd.simple.long %>% 
  
  # set up data
  filter(Model == 'recruitment probability') %>% 
  count(parameter, effect) %>% #get counts of + or - effects for each parameter
  mutate(parameter = factor(parameter)) %>% 
  mutate(effect = factor(effect)) %>% 
  mutate(perc = if_else(parameter %in% pp, n/2 * 100, #for params only used for 2 species
                        n/14 * 100)) %>% #for all other params used for 14 species
  
  # plot data
  ggplot(aes(fill=effect, y=perc, x=parameter)) +  
  geom_bar(position="dodge", stat="identity") +
  labs(x = NULL, y = NULL, tag = '(a)', title = paste('Recruitment (', n,'models )')) + 
  scale_fill_manual(values = cc2) + #add color-blind-friendly colour palette
  mytheme + theme(plot.tag.position = c(0.04, 0.985), legend.position = 'none') +
  ylim(0, 100)


# Recruitment counts

n <- dd.simple.long %>% 
  filter(Model == 'number of recruits' | Model == 'number of recruits*') %>% 
  distinct(Species) %>% 
  nrow() #number of models for title
## NOTE: only 13 models because best model for E. lanatum was intercept-only model

reccounts <- dd.simple.long %>% 
  filter(Model == 'number of recruits' | Model == 'number of recruits*') %>% #filter for model type
  count(parameter, effect) %>% #get counts of + or - effects for each parameter
  mutate(parameter = factor(parameter)) %>% #make parameter values factors
  mutate(effect = factor(effect)) %>% #make effect values factors
  mutate(perc = if_else(parameter %in% pp, n/1 * 100, #for params only used for 1 species (E. lanatum left out here, so only 1 species)
                        n/13 * 100)) %>% #for all other params used for 13 species
  ggplot(aes(fill=effect, y=perc, x=parameter)) +  
  geom_bar(position="dodge", stat="identity") +
  labs(x = NULL, y = NULL, tag = '(b)', title = paste('Recruit Counts (', n,'models )')) + 
  scale_fill_manual(values = cc2) + #add color-blind-friendly colour palette
  mytheme + theme(plot.tag.position = c(0.04, 0.985)) +
  ylim(0, 100)


# 1-year seedling survival

n <- dd.simple.long %>% 
  filter(Model == '1-yr survival') %>% 
  distinct(Species) %>% 
  nrow() #number of models for title

sdlg1yr <- dd.simple.long %>% 
  filter(Model == '1-yr survival') %>% #filter for model type
  count(parameter, effect) %>% #get counts of + or - effects for each parameter
  mutate(parameter = factor(parameter)) %>% #make parameter values factors
  mutate(effect = factor(effect)) %>% #make effect values factors
  mutate(perc = if_else(parameter %in% pp, n/1 * 100, #for params only used for 1 species
                        n/8 * 100)) %>% #for all other params used for 8 species
  ggplot(aes(fill=effect, y=perc, x=parameter)) +  
  geom_bar(position="dodge", stat="identity") +
  labs(x = NULL, y = NULL, tag = '(c)', title = paste('1-Year Survival (', n,'models )')) + 
  scale_fill_manual(values = cc2) + #add color-blind-friendly colour palette
  mytheme + theme(plot.tag.position = c(0.04, 0.985), legend.position = 'none') +
  ylim(0, 100)


# 2-year seedling survival

n <- dd.simple.long %>% 
  filter(Model == '2-yr survival*' | Model == '2-yr survival') %>% 
  distinct(Species) %>% 
  nrow() #number of models for title

sdlg2yr <- dd.simple.long %>% 
  filter(Model == '2-yr survival*' | Model == '2-yr survival') %>% #filter for model type
  count(parameter, effect) %>% #get counts of + or - effects for each parameter
  mutate(parameter = factor(parameter)) %>% #make parameter values factors
  mutate(effect = factor(effect)) %>% #make effect values factors
  mutate(perc = if_else(parameter %in% pp, n/1 * 100, #for params only used for 1 species
                        n/6 * 100)) %>% #for all other params used for 6 species
  ggplot(aes(fill=effect, y=perc, x=parameter)) +  
  geom_bar(position="dodge", stat="identity") +
  labs(x = NULL, y = NULL, tag = '(d)', title = paste('2-Year Survival (', n,'models )')) + 
  scale_fill_manual(values = cc2) + #add color-blind-friendly colour palette
  mytheme + theme(plot.tag.position = c(0.04, 0.985), legend.position = 'none') +
  ylim(0, 100)

# Combine plots and save
(fig <- grid.arrange(recprob, reccounts, sdlg1yr, sdlg2yr, nrow = 2,
                     left = textGrob("Proportion of Effect on Models [%]", #text for universal y-axis
                                     rot = 90, #rotate y-axis label 90 degrees
                                     vjust = 0.4, #adjust L-R position of universal y-axis
                                     gp = gpar(fontsize = 16)))) #set font size for universal y-axis

ggsave('outputs/ms_revisions_Dec2023/param_summaries.pdf', fig, width = 8, height = 6, units = 'in')




####################################################################################################

# # FIG. 4 [SUITABILITY-PREDS] # # 

####################################################################################################

# Data
load('outputs/ms_results_Jul2023/pred_bin.RData') #binary recruitment predictions at plot level (microhab_suitab.R)
load('outputs/ms_results_Jul2023/pred_cont.RData') #continuous recruitment predictions at plot level (microhab_suitab.R)
load('outputs/ms_results_Jul2023/pred_sdlg.RData') #seedling predictions at plot level (microhab_suitab.R)

# Theme
pp <- c(2, 5) #point size
tt <- 12 #text size
tt.leg <- 14 #legend & title size
cc <- c('#008000', 'black') #points color (https://html-color.codes/green)
ss <- c(16, 1) #hollow for beyond RL, solid for within # switched order to match with breaks

mytheme <- theme(axis.text.x = element_text(colour = 'black', size = tt),
                 axis.text.y = element_text(colour = 'black', size = tt),
                 plot.title = element_text(size = tt.leg), 
                 plot.tag = element_text(size = tt.leg), plot.tag.position = c(0.04, 0.98),
                 panel.background = element_blank(), axis.line = element_line(colour = "black"), 
                 legend.position = 'none'
) 

sp.list <- c('LUPLAT', 'TOLMEN', 'VACPAR') #manually chosen species to illustrate different suitab. patterns
spp <- pred.bin %>% #species names for each plot
  filter(species %in% sp.list) %>% 
  distinct(full_species) %>% 
  arrange(full_species) %>% 
  pull(full_species) #turn into vector


# RECRUITMENT PROB PLOTS
dat <- pred.bin %>% #generic dataframe 
  mutate(recruit = if_else(new_all == 0, 0, 1)) %>% #binary recruitment data 
  mutate(rl = if_else(temprange_pos < -1, 0, 1)) %>% #create binary thermal range limit column
  mutate(rl = as.factor(rl)) #needs to be factor for shape 

plot1 <- dat %>%
  filter(full_species == spp[1]) %>% 
  ggplot(aes(x = elev, y = pred.recr.bin, color = region)) + 
  geom_point(aes(size = recruit, shape = rl)) + 
  scale_shape_manual(values = ss, breaks = c(1,0)) +
  scale_size(range = pp) +
  scale_color_manual(values = cc) +
  labs(y = NULL, x = NULL, title = spp[1], tag = '(a)') +
  mytheme + theme(plot.title = element_text(face = 'italic'), legend.position = 'bottom') 

plot2 <- dat %>%
  filter(full_species == spp[2]) %>% 
  ggplot(aes(x = elev, y = pred.recr.bin, color = region)) + 
  geom_point(aes(size = recruit, shape = rl)) + 
  scale_shape_manual(values = ss, breaks = c(1,0)) +
  scale_size(range = pp) +
  scale_color_manual(values = cc) +
  labs(y = NULL, x = NULL, title = spp[2], tag = '(b)') +
  mytheme + theme(plot.title = element_text(face = 'italic'))

plot3 <- dat %>%
  filter(full_species == spp[3]) %>% 
  ggplot(aes(x = elev, y = pred.recr.bin, color = region)) + 
  geom_point(aes(size = recruit, shape = rl)) + 
  scale_shape_manual(values = ss, breaks = c(1,0)) +
  scale_size(range = pp) +
  scale_color_manual(values = cc) +
  labs(y = NULL, x = NULL, title = spp[3], tag = '(c)') +
  mytheme + theme(plot.title = element_text(face = 'italic'))


# RECRUITMENT COUNTS PLOTS
dat <- pred.cont %>% #generic dataframe 
  mutate(rl = if_else(temprange_pos < -1, 0, 1)) %>% #create binary thermal range limit column
  mutate(rl = as.factor(rl)) #needs to be factor for shape 

plot4 <- dat %>%
  filter(full_species == spp[1]) %>% 
  ggplot(aes(x = elev, y = pred.recr.cont, color = region)) + 
  geom_point(aes(size = rel_rec, shape = rl)) + 
  scale_shape_manual(values = ss, breaks = c(1,0)) +
  scale_size(range = pp) +
  scale_color_manual(values = cc) +
  labs(y = NULL, x = NULL, title = spp[1], tag = '(d)') +
  mytheme + theme(plot.title = element_text(face = 'italic'))

plot5 <- dat %>%
  filter(full_species == spp[2]) %>% 
  ggplot(aes(x = elev, y = pred.recr.cont, color = region)) + 
  geom_point(aes(size = rel_rec, shape = rl)) + 
  scale_shape_manual(values = ss, breaks = c(1,0)) +
  scale_size(range = pp) +
  scale_color_manual(values = cc) +
  labs(y = NULL, x = NULL, title = spp[2], tag = '(e)') +
  mytheme + theme(plot.title = element_text(face = 'italic'))

plot6 <- dat %>%
  filter(full_species == spp[3]) %>% 
  ggplot(aes(x = elev, y = pred.recr.cont, color = region)) + 
  geom_point(aes(size = rel_rec, shape = rl)) + 
  scale_shape_manual(values = ss, breaks = c(1,0)) +
  scale_size(range = pp) +
  scale_color_manual(values = cc) +
  labs(y = NULL, x = NULL, title = spp[3], tag = '(f)') +
  mytheme + theme(plot.title = element_text(face = 'italic'))


## 1-YEAR SEEDLING SURV
cc <- c('#008000', '#8fbc8f', 'black', '#a9a9a9') #points color (https://html-color.codes/green)

### KP added ###
pred.sdlg$region <- as.character(pred.sdlg$region)
####

dat <- pred.sdlg %>% 
  mutate(surv1yr = if_else(is.na(surv1yr), -1, surv1yr)) %>% #give numerical values to NA values for point size
  mutate(region = if_else(surv1yr < 0 & region == 'MB', 'MB-NA', region)) %>% #different region name for NA for color
  mutate(region = if_else(surv1yr < 0 & region == 'RP', 'RP-NA', region)) %>% 
  mutate(rl = if_else(temprange_pos < -1, 0, 1)) %>% #create binary thermal range limit column
  mutate(rl = as.factor(rl)) #needs to be factor for shape 

plot7 <- dat %>%
  filter(full_species == spp[1]) %>% 
  ggplot(aes(x = elev, y = pred.1yr.surv, color = region)) + 
  geom_point(aes(size = surv1yr, shape = rl)) + 
  scale_shape_manual(values = ss, breaks = c(1,0)) +
  scale_size(range = pp) +
  scale_color_manual(values = cc) +
  labs(y = NULL, x = NULL, title = spp[1], tag = '(g)') +
  mytheme + theme(plot.title = element_text(face = 'italic'))

plot8 <- dat %>%
  filter(full_species == spp[2]) %>% 
  ggplot(aes(x = elev, y = pred.1yr.surv, color = region)) + 
  geom_point(aes(size = surv1yr, shape = rl)) + 
  scale_shape_manual(values = ss, breaks = c(1,0)) +
  scale_size(range = pp) +
  scale_color_manual(values = cc) +
  labs(y = NULL, x = NULL, title = spp[2], tag = '(h)') +
  mytheme + theme(plot.title = element_text(face = 'italic'))

plot9 <- dat %>%
  filter(full_species == spp[3]) %>% 
  ggplot(aes(x = elev, y = pred.1yr.surv, color = region)) + 
  geom_point(aes(size = surv1yr, shape = rl)) + 
  scale_shape_manual(values = ss, breaks = c(1,0) , labels = c( "Below", "Beyond")) +
  scale_size(range = pp, guide = 'none') +
  scale_color_manual(values = cc, labels = c("West", "West-no recruitment", "East", "East-no recruitment")) +
  labs(y = NULL, x = NULL, title = spp[3], tag = '(i)') +
  mytheme + theme(plot.title = element_text(face = 'italic'), legend.position = 'bottom', legend.box="horizontal") +guides(color=guide_legend(title = 'Transect', title.position = 'top', nrow=1, byrow=TRUE,), shape = guide_legend(title = 'Range Limit Position', title.position = 'top', nrow = 1, byrow = TRUE)) 

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

plotLegend <- g_legend(plot9)
plot(plotLegend)
# Combine plots and save
fig <- grid.arrange(arrangeGrob(plot1+ theme(legend.position="none"), plot2, plot3, plot4, plot5, plot6, plot7 , plot8, plot9 + theme(legend.position = 'none'), heights=c(3,3, 3)), plotLegend, heights = c(10,1),
                    left = textGrob("Predicted Microhabitat Suitability", 
                                    rot = 90, gp = gpar(fontsize = tt.leg)), 
                    bottom = textGrob('Elevation [m]', gp = gpar(fontsize = tt.leg))
)

ggsave('outputs/ms_revisions_Dec2023//suitability-preds.pdf', fig, width = 8, height = 8, units = c('in'))




####################################################################################################

# # FIG. S1 [STUDY DESIGN] # # 

####################################################################################################

## Manually created on KG's computer.




####################################################################################################

# # FIG. S2 [SEED-SOURCE] # # 

####################################################################################################

# Note: germrate calculated by KG in scripts/range_position/prep/3_germ_rate_calc.R

# Data
load('data/tidy/census_spp.Rdata') # census data for species 2018-2022
seedDat <- read_csv('data/tidy/census_spp_seed_sources.csv')

# Change LUPLAT seed source to 'both' to correct data entry error from LUPLAT/LUPARC
seedDat <- seedDat %>% 
  mutate(seed_source = if_else(species == 'LUPLAT', 'both', seed_source))

# Combine DFs
census_plus_germ <- inner_join(census.spp, seedDat, by = "species")
census_plus_germ <- census_plus_germ %>% # reordering columns so seed_source is next to species
  relocate(seed_source, .after = species) %>%
  filter(germrate > 0)


# Theme
tt <- 16 #text size
cc <- c('grey', '#008000', 'darkblue') #points color (https://html-color.codes/green)

mytheme <- theme(axis.title.x = element_text(size = tt), axis.text.x = element_text(colour = 'black', size = tt),
                 axis.title.y = element_text(size = tt), axis.text.y = element_text(colour = 'black', size = tt),
                 panel.background = element_blank(), axis.line = element_line(colour = "black"),
                 plot.tag = element_text(size = tt), plot.tag.position = c(0.02, 0.98),
                 legend.text = element_text(size = tt), legend.title = element_blank(),
                 legend.key = element_rect(fill = "transparent"),
                 legend.position = c(0.65, 0.999), legend.justification = c(0, 1), 
                 legend.margin = margin(0, 6, 3, 0)
) 


# Histogram
germ_above_zero <- gghistogram(census_plus_germ, x = "germrate", y = "..density..", 
                               bins = 50, color = "seed_source", fill = "seed_source",
                               palette = c(cc[1], cc[2], cc[3])) +
  labs(x = 'Germination Rate', y = 'Frequency', tag = '(a)') +
  mytheme


# Boxplot
germination_boxplot <- ggplot(census_plus_germ, aes(seed_source, germrate, color=seed_source)) +
  geom_boxplot(position=position_dodge(1)) +
  labs(x ="Seed Source", y = "Germination Rate", tag = '(b)') +
  scale_color_manual(values = c(cc[1], cc[2], cc[3])) +
  mytheme + theme(legend.position = 'none')

# Combine plots and save
(fig <- grid.arrange(germ_above_zero, germination_boxplot, nrow = 1))

ggsave('outputs/ms_revisions_Dec2023/seed_source.pdf', fig, width = 8, height = 4, units = 'in')


# Histogram by species
tt <- 45

mytheme <- theme(axis.title.x = element_text(size = tt), axis.text.x = element_text(colour = 'black', size = tt),
                 axis.title.y = element_text(size = tt), axis.text.y = element_text(colour = 'black', size = tt),
                 panel.background = element_blank(), axis.line = element_line(colour = "black"),
                 plot.tag = element_text(size = tt), plot.tag.position = c(0.02, 0.98),
                 legend.text = element_text(size = tt), legend.title = element_blank(),
                 legend.key = element_rect(fill = "transparent"),
                 legend.position = c(0.65, 0.999), legend.justification = c(0, 1), 
                 legend.margin = margin(0, 6, 3, 0), plot.title = element_text(size=tt)
) 

v.spp <- unique(census_plus_germ$species) #vector of species
v.let <- letters[3:(length(v.spp)+2)] #vector for figure tags
v.sou <- c('both', 'local', 'nursery') #vector of seed sources

for (i in 1:length(v.spp)) {
  
  dat <- census_plus_germ %>% 
    filter(species == v.spp[i])
  
  gplot <- gghistogram(dat, x = "germrate", y = "..density..", 
                       bins = 50, color = "seed_source", fill = "seed_source",
                       palette = if(dat$seed_source[1] == 'both') {cc[1]} #color for seed source
                       else if(dat$seed_source[1] == 'local') {cc[2]}
                       else if(dat$seed_source[1] == 'nursery') {cc[3]}) +

    labs(x = 'Germination Rate', y = 'Frequency', tag = paste('(', v.let[i], ')'), 
         title = paste(v.spp[i])) +
    mytheme
  
  ggsave(paste0("outputs/ms_revisions_Dec2023/germrate", v.spp[i],".pdf"), gplot)
}


# Plots manually combined with (a) and (b) in PP on NC's computer

# Tukey-Kramer test to test which group means are different
TukeyHSD(aov(germrateANOVA))




####################################################################################################

# # FIG. S3 [PCA] # #

####################################################################################################

# Data 
load('data/tidy/census_spp.RData') #adds relative recruitment to census.tms (microhab_analysis_prep.R)

dat2 <- census.spp %>% 
  select(region, site, site1, species, rel_rec, canopycont, f_b, whc, c_n, tmoist_avgmin, #filter to env variables used in model selection
         t1_avgmax, winteravgmn, springtotalsnow, t3_avgmin) %>% 
  filter_at(vars(species, rel_rec, canopycont, f_b, whc, c_n, tmoist_avgmin,
                 t1_avgmax, winteravgmn, springtotalsnow, t3_avgmin), all_vars(!is.na(.)))

pp <- 4 #point size
tt <- 12 #text size
tt.leg <- 14 #legend & title size
cc <- c('#008000', 'black') #points color (https://html-color.codes/green)

mytheme <- theme(axis.title.x = element_text(size = tt), axis.text.x = element_text(size = tt, colour = 'black'),
                 axis.title.y = element_text(size = tt), axis.text.y = element_text(size = tt, colour = 'black'),
                 plot.title = element_text(size = tt.leg), 
                 plot.tag = element_text(size = tt.leg), plot.tag.position = c(0.04, 0.98),
                 panel.background = element_blank(), axis.line = element_line(colour = "black"),
                 legend.position = "none"
) 

## Entire microhabitat space per plot ##
# source: https://www.statology.org/principal-components-analysis-in-r/

# Calculate principal components
results2 <- prcomp(dat2[, -c(1:5)], scale = TRUE)

# Reverse the signs
results2$rotation <- -1*results2$rotation

# Display principal components
results2$rotation

# Reverse the signs of the scores
results2$x <- -1*results2$x

# Display the first six scores
head(results2$x)

# # Calculate proportion of variance explained by PC1 and PC2
# ss <- c(summary(results)$importance[2, 1], summary(results)$importance[2, 2]) #prop var for PC1, PC2
# xx <- paste('PC1 (', round(ss[1], 2) * 100, '% )') #axes labels
# yy <- paste('PC2 (', round(ss[2], 2) * 100, '% )')
# summary(results)$importance #prop var for all axes => 84% prop variances with PC1-5


# Plot results

dat2$site2 <- do.call(rbind, strsplit(as.character(dat2$site1), split = "-"))[,2]
dat2$site2[which(dat2$site2 == "1B")] <- "1"

pca_plot <- autoplot(results2, loadings = T,  loadings.label = F, loadings.label.size = 6, data = dat2, 
         colour = 'region', loadings.colour = "black", size = 3, loadings.label.color = 'black', 
         loadings.label.vjust=-0.5, loadings.label.hjust = 0.1) + 
  mytheme + scale_shape_manual(values = c(1,2)) + 
  scale_color_manual(values = cc) # with color = region

# autoplot(results2, loadings = T,  loadings.label = TRUE, loadings.label.size = 6, data = dat2, 
#          colour = 'site2', shape = 'region', loadings.colour = "black", size = 3, 
#          loadings.label.color = 'black', loadings.label.vjust=-0.5, loadings.label.hjust = 0.1) + 
#   mytheme + scale_shape_manual(values = c(1,2))  #with shape = region and color = 1-15 site in each transect


pdf('outputs/ms_revisions_Apr2024/allspeciesPCA_colorRegion.pdf')
pca_plot
dev.off()




####################################################################################################

# # FIG. S4 [PARAM-SPP] # # 

####################################################################################################

# Data
load('outputs/ms_results_Jul2023/ALL_mod_avg_SIMPLE.RData') #model averaging results for all models (ms_tables.R)

# Number of times each parameter used in model in long DF format 
dd.simple.long <- gather(dd.simple, parameter, effect, 'Year':'Spring Days with Snow^2', factor_key=TRUE) #making dd.tab.simple df long so that parameter values can be grouped by +/- in stacked bar graph
is.na(dd.simple.long$effect) <- dd.simple.long$effect == 0 #make 0 values into NAs
dd.simple.long <- drop_na(dd.simple.long) #drop NA values from df

# Change parameter names to be shorter in figures
params <- unique(dd.simple.long$parameter)
dd.simple.long <- dd.simple.long %>%
  mutate(parameter = if_else(parameter == as.character(params[3]), 'Canopy', parameter)) %>%
  mutate(parameter = if_else(parameter == as.character(params[4]), 'Canopy^2', parameter)) %>%
  mutate(parameter = if_else(parameter == as.character(params[5]), 'F:B', parameter)) %>%
  mutate(parameter = if_else(parameter == as.character(params[6]), 'F:B^2', parameter)) %>%
  mutate(parameter = if_else(parameter == as.character(params[7]), 'C:N', parameter)) %>%
  mutate(parameter = if_else(parameter == as.character(params[8]), 'C:N^2', parameter)) %>%
  mutate(parameter = if_else(parameter == as.character(params[9]), 'WHC', parameter)) %>%
  mutate(parameter = if_else(parameter == as.character(params[10]), 'WHC^2', parameter)) %>%
  mutate(parameter = if_else(parameter == as.character(params[11]), 'S Temp', parameter)) %>%
  mutate(parameter = if_else(parameter == as.character(params[12]), 'S Soil Moist', parameter)) %>%
  mutate(parameter = if_else(parameter == as.character(params[13]), 'W Soil Temp', parameter)) %>%
  mutate(parameter = if_else(parameter == as.character(params[14]), 'W Soil Temp^2', parameter)) %>%
  mutate(parameter = if_else(parameter == as.character(params[15]), 'Snow', parameter)) %>%
  mutate(parameter = if_else(parameter == as.character(params[16]), 'Snow^2', parameter))

# Vectors of species names and number of models run per species
spp <- as.data.frame(dd.simple %>% count(Species) %>% mutate(n = as.numeric(n)))

# Vector of TOMST parameters for E. lanatum (all mods), M. nervosa (all mods), R. ursinus (1/2 mods)
# Note: don't need to manually change anything for E. lan and M. nerv since TOMST params used in all models
pp <- c('S Temp', 'S Soil Moist')

# Theme
tt <- 15
tt.x <- 12
mytheme <- theme(axis.text.x = element_text(colour = 'black', size = tt.x, angle = 60, margin = margin(t = 10), hjust = 1),
                 axis.title.y = element_text(size = tt), axis.text.y = element_text(colour = 'black', size = tt),
                 plot.title = element_text(size = tt), 
                 plot.tag = element_text(size = tt), plot.tag.position = c(0.1, 0.97),
                 panel.background = element_blank(), axis.line = element_line(colour = "black"),
                 legend.position = 'none', legend.margin = margin(0, 6, 3, 0), 
                 legend.box.background = element_rect(colour = 'black'))

cc2 <- c('#AA4499', '#008000') #color-blind friendly palette to use for +/- effect 


# Species plots
one <- dd.simple.long %>% 
  filter(Species == spp$Species[1]) %>% #***
  count(parameter, effect) %>% 
  mutate(parameter = factor(parameter)) %>% #make parameter values factors
  mutate(effect = factor(effect)) %>% #make effect values factors
  mutate(perc = n/spp$n[1] * 100) %>%  #n/x with x = total models run for the species #***
  ggplot(aes(fill=effect, y=perc, x=parameter)) +  
  geom_bar(position="dodge", stat="identity") +
  labs(x = NULL, y = NULL, tag = '(a)', title = paste(spp$Species[1])) + #***
  scale_fill_manual(values = cc2) + #add color-blind-friendly colour palette
  mytheme + theme(plot.tag.position = c(0.04, 0.985)) +
  mytheme + theme(plot.title = element_text(face = 'italic')) +
  theme(legend.text = element_text(size = tt.x), legend.title = element_blank(),
        legend.key = element_rect(fill = "transparent"),
        legend.position = c(0.9, 1)) +
  ylim(0, 100)

two <- dd.simple.long %>% 
  filter(Species == spp$Species[2]) %>% 
  count(parameter, effect) %>% 
  mutate(parameter = factor(parameter)) %>% #make parameter values factors
  mutate(effect = factor(effect)) %>% #make effect values factors
  mutate(perc = n/spp$n[2] * 100) %>%  #n/x with x = total models run for the species 
  ggplot(aes(fill=effect, y=perc, x=parameter)) +  
  geom_bar(position="dodge", stat="identity") +
  labs(x = NULL, y = NULL, tag = '(b)', title = paste(spp$Species[2])) + 
  scale_fill_manual(values = cc2) + #add color-blind-friendly colour palette
  mytheme + theme(plot.tag.position = c(0.04, 0.985)) +
  mytheme + theme(plot.title = element_text(face = 'italic')) +
  ylim(0, 100)

three <- dd.simple.long %>% 
  filter(Species == spp$Species[3]) %>% 
  count(parameter, effect) %>% 
  mutate(parameter = factor(parameter)) %>% #make parameter values factors
  mutate(effect = factor(effect)) %>% #make effect values factors
  mutate(perc = n/spp$n[3] * 100) %>%  #n/x with x = total models run for the species 
  ggplot(aes(fill=effect, y=perc, x=parameter)) +  
  geom_bar(position="dodge", stat="identity") +
  labs(x = NULL, y = NULL, tag = '(c)', title = paste(spp$Species[3])) + 
  scale_fill_manual(values = cc2) + #add color-blind-friendly colour palette
  mytheme + theme(plot.tag.position = c(0.04, 0.985)) +
  mytheme + theme(plot.title = element_text(face = 'italic')) +
  ylim(0, 100)

four <- dd.simple.long %>% 
  filter(Species == spp$Species[4]) %>% 
  count(parameter, effect) %>% 
  mutate(parameter = factor(parameter)) %>% #make parameter values factors
  mutate(effect = factor(effect)) %>% #make effect values factors
  mutate(perc = n/spp$n[4] * 100) %>%  #n/x with x = total models run for the species 
  ggplot(aes(fill=effect, y=perc, x=parameter)) +  
  geom_bar(position="dodge", stat="identity") +
  labs(x = NULL, y = NULL, tag = '(d)', title = paste(spp$Species[4])) + 
  scale_fill_manual(values = cc2) + #add color-blind-friendly colour palette
  mytheme + theme(plot.tag.position = c(0.04, 0.985)) +
  mytheme + theme(plot.title = element_text(face = 'italic')) +
  ylim(0, 100)

five <- dd.simple.long %>% 
  filter(Species == spp$Species[5]) %>% 
  count(parameter, effect) %>% 
  mutate(parameter = factor(parameter)) %>% #make parameter values factors
  mutate(effect = factor(effect)) %>% #make effect values factors
  mutate(perc = n/spp$n[5] * 100) %>%  #n/x with x = total models run for the species 
  ggplot(aes(fill=effect, y=perc, x=parameter)) +  
  geom_bar(position="dodge", stat="identity") +
  labs(x = NULL, y = NULL, tag = '(e)', title = paste(spp$Species[5])) + 
  scale_fill_manual(values = cc2) + #add color-blind-friendly colour palette
  mytheme + theme(plot.tag.position = c(0.04, 0.985)) +
  mytheme + theme(plot.title = element_text(face = 'italic')) +
  ylim(0, 100)

six <- dd.simple.long %>% 
  filter(Species == spp$Species[6]) %>% 
  count(parameter, effect) %>% 
  mutate(parameter = factor(parameter)) %>% #make parameter values factors
  mutate(effect = factor(effect)) %>% #make effect values factors
  mutate(perc = n/spp$n[6] * 100) %>%  #n/x with x = total models run for the species 
  ggplot(aes(fill=effect, y=perc, x=parameter)) +  
  geom_bar(position="dodge", stat="identity") +
  labs(x = NULL, y = NULL, tag = '(f)', title = paste(spp$Species[6])) + 
  scale_fill_manual(values = cc2) + #add color-blind-friendly colour palette
  mytheme + theme(plot.tag.position = c(0.04, 0.985)) +
  mytheme + theme(plot.title = element_text(face = 'italic')) +
  ylim(0, 100)

seven <- dd.simple.long %>% 
  filter(Species == spp$Species[7]) %>% 
  count(parameter, effect) %>% 
  mutate(parameter = factor(parameter)) %>% #make parameter values factors
  mutate(effect = factor(effect)) %>% #make effect values factors
  mutate(perc = n/spp$n[7] * 100) %>% 
  ggplot(aes(fill=effect, y=perc, x=parameter)) +  
  geom_bar(position="dodge", stat="identity") +
  labs(x = NULL, y = NULL, tag = '(g)', title = paste(spp$Species[7])) + 
  scale_fill_manual(values = cc2) + #add color-blind-friendly colour palette
  mytheme + theme(plot.tag.position = c(0.04, 0.985)) +
  mytheme + theme(plot.title = element_text(face = 'italic')) +
  ylim(0, 100)

eight <- dd.simple.long %>% 
  filter(Species == spp$Species[8]) %>% 
  count(parameter, effect) %>% 
  mutate(parameter = factor(parameter)) %>% #make parameter values factors
  mutate(effect = factor(effect)) %>% #make effect values factors
  mutate(perc = n/spp$n[8] * 100) %>%  #n/x with x = total models run for the species 
  ggplot(aes(fill=effect, y=perc, x=parameter)) +  
  geom_bar(position="dodge", stat="identity") +
  labs(x = NULL, y = NULL, tag = '(h)', title = paste(spp$Species[8])) + 
  scale_fill_manual(values = cc2) + #add color-blind-friendly colour palette
  mytheme + theme(plot.tag.position = c(0.04, 0.985)) +
  mytheme + theme(plot.title = element_text(face = 'italic')) +
  ylim(0, 100)

# Manual percentage corrections for R ursinus
nine <- dd.simple.long %>% 
  filter(Species == spp$Species[9]) %>% 
  count(parameter, effect) %>% 
  mutate(parameter = factor(parameter)) %>% #make parameter values factors
  mutate(effect = factor(effect)) %>% #make effect values factors
  mutate(perc = if_else(parameter %in% pp, n/spp$n[9] *2 * 100, #for params only used for half of R ursinus models
                        n/spp$n[9] * 100)) %>%  #n/x with x = total models run for the species 
  ggplot(aes(fill=effect, y=perc, x=parameter)) +  
  geom_bar(position="dodge", stat="identity") +
  labs(x = NULL, y = NULL, tag = '(i)', title = paste(spp$Species[9])) + 
  scale_fill_manual(values = cc2) + #add color-blind-friendly colour palette
  mytheme + theme(plot.tag.position = c(0.04, 0.985)) +
  mytheme + theme(plot.title = element_text(face = 'italic')) +
  ylim(0, 100)

ten <- dd.simple.long %>% 
  filter(Species == spp$Species[10]) %>% 
  count(parameter, effect) %>% 
  mutate(parameter = factor(parameter)) %>% #make parameter values factors
  mutate(effect = factor(effect)) %>% #make effect values factors
  mutate(perc = n/spp$n[10] * 100) %>%  #for all other params used in all models
  ggplot(aes(fill=effect, y=perc, x=parameter)) +  
  geom_bar(position="dodge", stat="identity") +
  labs(x = NULL, y = NULL, tag = '(j)', title = paste(spp$Species[10])) + 
  scale_fill_manual(values = cc2) + #add color-blind-friendly colour palette
  mytheme + theme(plot.tag.position = c(0.04, 0.985)) +
  mytheme + theme(plot.title = element_text(face = 'italic')) +
  ylim(0, 100)

eleven <- dd.simple.long %>% 
  filter(Species == spp$Species[11]) %>% 
  count(parameter, effect) %>% 
  mutate(parameter = factor(parameter)) %>% #make parameter values factors
  mutate(effect = factor(effect)) %>% #make effect values factors
  mutate(perc = n/spp$n[11] * 100) %>%  #n/x with x = total models run for the species 
  ggplot(aes(fill=effect, y=perc, x=parameter)) +  
  geom_bar(position="dodge", stat="identity") +
  labs(x = NULL, y = NULL, tag = '(k)', title = paste(spp$Species[11])) + 
  scale_fill_manual(values = cc2) + #add color-blind-friendly colour palette
  mytheme + theme(plot.tag.position = c(0.04, 0.985)) +
  mytheme + theme(plot.title = element_text(face = 'italic')) +
  ylim(0, 100)

twelve <- dd.simple.long %>% 
  filter(Species == spp$Species[12]) %>% 
  count(parameter, effect) %>% 
  mutate(parameter = factor(parameter)) %>% #make parameter values factors
  mutate(effect = factor(effect)) %>% #make effect values factors
  mutate(perc = n/spp$n[12] * 100) %>%  #for all other params used in all models
  ggplot(aes(fill=effect, y=perc, x=parameter)) +  
  geom_bar(position="dodge", stat="identity") +
  labs(x = NULL, y = NULL, tag = '(l)', title = paste(spp$Species[12])) + 
  scale_fill_manual(values = cc2) + #add color-blind-friendly colour palette
  mytheme + theme(plot.tag.position = c(0.04, 0.985)) +
  mytheme + theme(plot.title = element_text(face = 'italic')) +
  ylim(0, 100)

thirt <- dd.simple.long %>% 
  filter(Species == spp$Species[13]) %>% 
  count(parameter, effect) %>% 
  mutate(parameter = factor(parameter)) %>% #make parameter values factors
  mutate(effect = factor(effect)) %>% #make effect values factors
  mutate(perc = n/spp$n[13] * 100) %>%  #n/x with x = total models run for the species 
  ggplot(aes(fill=effect, y=perc, x=parameter)) +  
  geom_bar(position="dodge", stat="identity") +
  labs(x = NULL, y = NULL, tag = '(m)', title = paste(spp$Species[13])) + 
  scale_fill_manual(values = cc2) + #add color-blind-friendly colour palette
  mytheme + theme(plot.tag.position = c(0.04, 0.985)) +
  mytheme + theme(plot.title = element_text(face = 'italic')) +
  ylim(0, 100)

fourt <- dd.simple.long %>% 
  filter(Species == spp$Species[14]) %>% 
  count(parameter, effect) %>% 
  mutate(parameter = factor(parameter)) %>% #make parameter values factors
  mutate(effect = factor(effect)) %>% #make effect values factors
  mutate(perc = n/spp$n[14] * 100) %>%  #n/x with x = total models run for the species 
  ggplot(aes(fill=effect, y=perc, x=parameter)) +  
  geom_bar(position="dodge", stat="identity") +
  labs(x = NULL, y = NULL, tag = '(n)', title = paste(spp$Species[14])) + 
  scale_fill_manual(values = cc2) + #add color-blind-friendly colour palette
  mytheme + theme(plot.tag.position = c(0.04, 0.985)) +
  mytheme + theme(plot.title = element_text(face = 'italic')) +
  ylim(0, 100)

# Combine plots and save
fig <- grid.arrange(one, two, three, four, five, six, seven, eight, nine, ten, eleven, twelve, thirt, 
                    fourt, nrow = 5,
                    left = textGrob("Proportion of Effect on Models [%]", #text for universal y-axis
                                    rot = 90, #rotate y-axis label 90 degrees
                                    vjust = 0.4, #adjust L-R position of universal y-axis
                                    gp = gpar(fontsize = 16))) #set font size for universal y-axis


ggsave('outputs/ms_revisions_Dec2023/param_spp.pdf', fig, width = 8, height = 12, units = 'in')




####################################################################################################

# # FIG. S5 [SUITABILITY-RECRUITMENTPROB] # # 

####################################################################################################

# Data
load('outputs/ms_results_Jul2023/pred_bin.RData') #binary recruitment predictions at plot level (microhab_suitab.R)

# Theme
pp <- c(1, 3) #point size
tt <- 10 #text size
cc <- c('#008000', 'black') #points color (https://html-color.codes/green)
ss <- c(16, 1) #hollow for beyond RL, solid for within

mytheme <- theme(axis.title.x = element_text(size = tt), axis.text.x = element_text(colour = 'black', size = tt),
                 axis.title.y = element_text(size = tt), axis.text.y = element_text(colour = 'black', size = tt),
                 plot.title = element_text(size = tt), 
                 plot.tag = element_text(size = tt), plot.tag.position = c(0.05, 0.97),
                 panel.background = element_blank(), axis.line = element_line(colour = "black"),
                 legend.position = 'none'
) 


# RECRUITMENT PROB PLOTS
dat <- pred.bin %>% #generic dataframe 
  mutate(recruit = if_else(new_all == 0, 0, 1)) %>% #binary recruitment data 
  mutate(rl = if_else(temprange_pos < -1, 0, 1)) %>% #create binary thermal range limit column
  mutate(rl = as.factor(rl)) #needs to be factor for shape 

spp <- sort(unique(dat$full_species)) #species names for each plot

plot1 <- dat %>%
  filter(full_species == spp[1]) %>% 
  ggplot(aes(x = elev, y = pred.recr.bin, color = region)) + 
  geom_point(aes(size = recruit, shape = rl)) + 
  scale_shape_manual(values = ss, breaks= c(1,0)) +
  scale_size(range = pp) +
  scale_color_manual(values = cc) +
  labs(y = NULL, x = NULL, title = spp[1], tag = '(a)') +
  mytheme + 
  theme(plot.title = element_text(face = 'italic')) 
  
plot2 <- dat %>%
  filter(full_species == spp[2]) %>% 
  ggplot(aes(x = elev, y = pred.recr.bin, color = region)) + 
  geom_point(aes(size = recruit, shape = rl)) + 
  scale_shape_manual(values = ss, breaks= c(1,0)) +
  scale_size(range = pp) +
  scale_color_manual(values = cc) +
  labs(y = NULL, x = NULL, title = spp[2], tag = '(b)') +
  mytheme + theme(plot.title = element_text(face = 'italic'))

plot3 <- dat %>%
  filter(full_species == spp[3]) %>% 
  ggplot(aes(x = elev, y = pred.recr.bin, color = region)) + 
  geom_point(aes(size = recruit, shape = rl)) + 
  scale_shape_manual(values = ss, breaks= c(1,0)) +
  scale_size(range = pp) +
  scale_color_manual(values = cc) +
  labs(y = NULL, x = NULL, title = spp[3], tag = '(c)') +
  mytheme + theme(plot.title = element_text(face = 'italic'))

plot4 <- dat %>%
  filter(full_species == spp[4]) %>% 
  ggplot(aes(x = elev, y = pred.recr.bin, color = region)) + 
  geom_point(aes(size = recruit, shape = rl)) + 
  scale_shape_manual(values = ss, breaks= c(1,0)) +
  scale_size(range = pp) +
  scale_color_manual(values = cc) +
  labs(y = NULL, x = NULL, title = spp[4], tag = '(d)') +
  mytheme + theme(plot.title = element_text(face = 'italic'))

plot5 <- dat %>%
  filter(full_species == spp[5]) %>% 
  ggplot(aes(x = elev, y = pred.recr.bin, color = region)) + 
  geom_point(aes(size = recruit, shape = rl)) + 
  scale_shape_manual(values = ss, breaks= c(1,0)) +
  scale_size(range = pp) +
  scale_color_manual(values = cc) +
  labs(y = NULL, x = NULL, title = spp[5], tag = '(e)') +
  mytheme + theme(plot.title = element_text(face = 'italic'))

plot6 <- dat %>%
  filter(full_species == spp[6]) %>% 
  ggplot(aes(x = elev, y = pred.recr.bin, color = region)) + 
  geom_point(aes(size = recruit, shape = rl)) + 
  scale_shape_manual(values = ss, breaks= c(1,0)) +
  scale_size(range = pp) +
  scale_color_manual(values = cc) +
  labs(y = NULL, x = NULL, title = spp[6], tag = '(f)') +
  mytheme + theme(plot.title = element_text(face = 'italic'))

plot7 <- dat %>%
  filter(full_species == spp[7]) %>% 
  ggplot(aes(x = elev, y = pred.recr.bin, color = region)) + 
  geom_point(aes(size = recruit, shape = rl)) + 
  scale_shape_manual(values = ss, breaks= c(1,0)) +
  scale_size(range = pp) +
  scale_color_manual(values = cc) +
  labs(y = NULL, x = NULL, title = spp[7], tag = '(g)') +
  mytheme + theme(plot.title = element_text(face = 'italic'))
 
plot8 <- dat %>%
  filter(full_species == spp[8]) %>% 
  ggplot(aes(x = elev, y = pred.recr.bin, color = region)) + 
  geom_point(aes(size = recruit, shape = rl)) + 
  scale_shape_manual(values = ss, breaks= c(1,0), labels = c("Below", "Beyond")) +
  scale_size(range = pp, guide = 'none') +
  scale_color_manual(values = cc, labels = c("West", "East")) +
  labs(y = NULL, x = NULL, title = spp[8], tag = '(h)') +
  mytheme + theme(plot.title = element_text(face = 'italic'), legend.position = 'bottom') +guides(color=guide_legend(title = 'Transect', title.position = 'top', nrow=1, byrow=TRUE,), shape = guide_legend(title = 'Range Limit Position', title.position = 'top', nrow = 1, byrow = TRUE)) 

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

plotLegend <- g_legend(plot8)
#plot(plotLegend)

# Combine plots and save
fig <- grid.arrange(arrangeGrob(plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8 + theme(legend.position = 'none'), heights=c(3,3)), plotLegend, heights = c(12,2),
                    left = textGrob("Predicted Microhabitat Suitability\n (recruitment)", #text for universal y-axis
                                    rot = 90, #rotate y-axis label 90 degrees
                                    gp = gpar(fontsize = tt)), #set font size for universal y-axis
                    bottom = textGrob('Elevation [m]', gp = gpar(fontsize = tt))
)


ggsave('outputs/ms_revisions_Dec2023/suitability-recruitprob.pdf', fig, width = 8, height = 4.75, units = c('in'))




####################################################################################################

# # FIG. S6 [SUITABILITY-RECRUITMENTCOUNT] # # 

####################################################################################################

# Data
load('outputs/ms_results_Jul2023/pred_cont.RData') #continuous recruitment predictions at plot level (microhab_suitab.R)

# Theme
pp <- c(1, 3) #point size
tt <- 10 #text size
cc <- c('#008000', 'black') #points color (https://html-color.codes/green)
ss <- c(1, 16) #hollow for beyond RL, solid for within

mytheme <- theme(axis.title.x = element_text(size = tt), axis.text.x = element_text(colour = 'black', size = tt),
                 axis.title.y = element_text(size = tt), axis.text.y = element_text(colour = 'black', size = tt),
                 plot.title = element_text(size = tt), 
                 plot.tag = element_text(size = tt), plot.tag.position = c(0.05, 0.97),
                 panel.background = element_blank(), axis.line = element_line(colour = "black"),
                 legend.position = 'none'
) 


# RECRUITMENT COUNTS PLOTS
dat <- pred.cont %>% #generic dataframe 
  mutate(rl = if_else(temprange_pos < -1, 0, 1)) %>% #create binary thermal range limit column
  mutate(rl = as.factor(rl)) #needs to be factor for shape 
spp <- sort(unique(dat$full_species)) #species names for each plot

plot1 <- dat %>%
  filter(full_species == spp[1]) %>% 
  ggplot(aes(x = elev, y = pred.recr.cont, color = region)) + 
  geom_point(aes(size = rel_rec, shape = rl)) + 
  scale_shape_manual(values = ss) +
  scale_size(range = pp) +
  scale_color_manual(values = cc) +
  labs(y = NULL, x = NULL, title = spp[1], tag = '(a)') +
  mytheme + theme(plot.title = element_text(face = 'italic'))

plot2 <- dat %>%
  filter(full_species == spp[2]) %>% 
  ggplot(aes(x = elev, y = pred.recr.cont, color = region)) + 
  geom_point(aes(size = rel_rec, shape = rl)) + 
  scale_shape_manual(values = ss) +
  scale_size(range = pp) +
  scale_color_manual(values = cc) +
  labs(y = NULL, x = NULL, title = spp[2], tag = '(b)') +
  mytheme + theme(plot.title = element_text(face = 'italic'))

plot3 <- dat %>%
  filter(full_species == spp[3]) %>% 
  ggplot(aes(x = elev, y = pred.recr.cont, color = region)) + 
  geom_point(aes(size = rel_rec, shape = rl)) + 
  scale_shape_manual(values = ss) +
  scale_size(range = pp) +
  scale_color_manual(values = cc) +
  labs(y = NULL, x = NULL, title = spp[3], tag = '(c)') +
  mytheme + theme(plot.title = element_text(face = 'italic'))

plot4 <- dat %>%
  filter(full_species == spp[4]) %>% 
  ggplot(aes(x = elev, y = pred.recr.cont, color = region)) + 
  geom_point(aes(size = rel_rec, shape = rl)) + 
  scale_shape_manual(values = ss) +
  scale_size(range = pp) +
  scale_color_manual(values = cc) +
  labs(y = NULL, x = NULL, title = spp[4], tag = '(d)') +
  mytheme + theme(plot.title = element_text(face = 'italic'))

plot5 <- dat %>%
  filter(full_species == spp[5]) %>% 
  ggplot(aes(x = elev, y = pred.recr.cont, color = region)) + 
  geom_point(aes(size = rel_rec, shape = rl)) + 
  scale_shape_manual(values = ss) +
  scale_size(range = pp) +
  scale_color_manual(values = cc) +
  labs(y = NULL, x = NULL, title = spp[5], tag = '(e)') +
  mytheme + theme(plot.title = element_text(face = 'italic'))

plot6 <- dat %>%
  filter(full_species == spp[6]) %>% 
  ggplot(aes(x = elev, y = pred.recr.cont, color = region)) + 
  geom_point(aes(size = rel_rec, shape = rl)) + 
  scale_shape_manual(values = ss) +
  scale_size(range = pp) +
  scale_color_manual(values = cc) +
  labs(y = NULL, x = NULL, title = spp[6], tag = '(f)') +
  mytheme + theme(plot.title = element_text(face = 'italic'))

plot7 <- dat %>%
  filter(full_species == spp[7]) %>% 
  ggplot(aes(x = elev, y = pred.recr.cont, color = region)) + 
  geom_point(aes(size = rel_rec, shape = rl)) + 
  scale_shape_manual(values = ss) +
  scale_size(range = pp) +
  scale_color_manual(values = cc) +
  labs(y = NULL, x = NULL, title = spp[7], tag = '(g)') +
  mytheme + theme(plot.title = element_text(face = 'italic'))


plot8 <- dat %>%
  filter(full_species == spp[8]) %>% 
  ggplot(aes(x = elev, y = pred.recr.cont, color = region)) + 
  geom_point(aes(size = rel_rec, shape = rl)) + 
  scale_shape_manual(values = ss, labels = c("Beyond", "Below")) +
  scale_size(range = pp, guide = 'none') +
  scale_color_manual(values = cc, labels = c("West", "East")) +
  labs(y = NULL, x = NULL, title = spp[8], tag = '(h)') +
  mytheme + theme(plot.title = element_text(face = 'italic'), legend.position = 'bottom') + guides(color=guide_legend(title = 'Transect', title.position = 'top', nrow=1, byrow=TRUE,), shape = guide_legend(title = 'Range Limit Position', title.position = 'top', nrow = 1, byrow = TRUE)) 

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

plotLegend <- g_legend(plot8)
plot(plotLegend)

fig <- grid.arrange(arrangeGrob(plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8 + theme(legend.position = 'none'), heights=c(3,3)), plotLegend, heights = c(10,1),
                    left = textGrob("Predicted Microhabitat Suitability\n (recruitment)", #text for universal y-axis
                                    rot = 90, #rotate y-axis label 90 degrees
                                    gp = gpar(fontsize = tt)), #set font size for universal y-axis
                    bottom = textGrob('Elevation [m]', gp = gpar(fontsize = tt))
)

# Combine plots and save
fig <- grid.arrange(arrangeGrob(plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8 + theme(legend.position = 'none'), heights=c(3,3)), plotLegend, heights = c(10,1), 
                     left = textGrob("Predicted Microhabitat Suitability\n (relative recruitment counts)", #text for universal y-axis
                                     rot = 90, #rotate y-axis label 90 degrees
                                     gp = gpar(fontsize = tt)), #set font size for universal y-axis
                     bottom = textGrob('Elevation [m]', gp = gpar(fontsize = tt))
)

ggsave('outputs/ms_revisions_Dec2023/suitability-recruitcounts.pdf', fig,  width = 8, height = 4.75, units = c('in'))




####################################################################################################

# # FIG. S7 [SUITABILITY-SDLGSURV] # # 

####################################################################################################

# Data
load('outputs/ms_results_Jul2023/pred_sdlg.RData') #seedling predictions at plot level (microhab_suitab.R)

# Theme
pp <- c(1, 3) #point size
tt <- 12 #text size
ss <- c(16,1) #hollow for beyond RL, solid for within

mytheme <- theme(axis.title.x = element_text(size = tt), axis.text.x = element_text(colour = 'black', size = tt),
                 axis.title.y = element_text(size = tt), axis.text.y = element_text(colour = 'black', size = tt),
                 plot.title = element_text(size = tt), 
                 plot.tag = element_text(size = tt), plot.tag.position = c(0.05, 0.97),
                 panel.background = element_blank(), axis.line = element_line(colour = "black"),
                 legend.position = 'none'
) 


## 1-YEAR SEEDLING SURV
cc <- c('#008000', '#8fbc8f', 'black', '#a9a9a9') #points color (https://html-color.codes/green)
pred.sdlg$region <- as.character(pred.sdlg$region)
dat <- pred.sdlg %>% 
  mutate(surv1yr = if_else(is.na(surv1yr), -1, surv1yr)) %>% #give numerical values to NA values for point size
  mutate(region = if_else(surv1yr < 0 & region == 'MB', 'MB-NA', region)) %>% #different region name for NA for color
  mutate(region = if_else(surv1yr < 0 & region == 'RP', 'RP-NA', region)) %>% 
  mutate(rl = if_else(temprange_pos < -1, 0, 1)) %>% #create binary thermal range limit column
  mutate(rl = as.factor(rl)) #needs to be factor for shape 

spp <- sort(unique(dat$full_species)) #species names for each plot

plot1 <- dat %>%
  filter(full_species == spp[1]) %>% 
  ggplot(aes(x = elev, y = pred.1yr.surv, color = region)) + 
  geom_point(aes(size = surv1yr, shape = rl)) + 
  scale_shape_manual(values = ss, breaks = c(1,0)) +
  scale_size(range = pp) +
  scale_color_manual(values = cc) +
  labs(y = NULL, x = NULL, title = spp[1], tag = '(a)') +
  mytheme + theme(plot.title = element_text(face = 'italic'))

plot2 <- dat %>%
  filter(full_species == spp[2]) %>% 
  ggplot(aes(x = elev, y = pred.1yr.surv, color = region)) + 
  geom_point(aes(size = surv1yr, shape = rl)) + 
  scale_shape_manual(values = ss, breaks = c(1,0)) +
  scale_size(range = pp) +
  scale_color_manual(values = cc) +
  labs(y = NULL, x = NULL, title = spp[2], tag = '(b)') +
  mytheme + theme(plot.title = element_text(face = 'italic'))

plot3 <- dat %>%
  filter(full_species == spp[3]) %>% 
  ggplot(aes(x = elev, y = pred.1yr.surv, color = region)) + 
  geom_point(aes(size = surv1yr, shape = rl)) + 
  scale_shape_manual(values = ss, breaks = c(1,0)) +
  scale_size(range = pp) +
  scale_color_manual(values = cc) +
  labs(y = NULL, x = NULL, title = spp[3], tag = '(c)') +
  mytheme + theme(plot.title = element_text(face = 'italic'))

plot4 <- dat %>%
  filter(full_species == spp[4]) %>% 
  ggplot(aes(x = elev, y = pred.1yr.surv, color = region)) + 
  geom_point(aes(size = surv1yr, shape = rl)) + 
  scale_shape_manual(values = ss, breaks = c(1,0)) +
  scale_size(range = pp) +
  scale_color_manual(values = cc) +
  labs(y = NULL, x = NULL, title = spp[4], tag = '(d)') +
  mytheme + theme(plot.title = element_text(face = 'italic'))

plot5 <- dat %>%
  filter(full_species == spp[5]) %>% 
  ggplot(aes(x = elev, y = pred.1yr.surv, color = region)) + 
  geom_point(aes(size = surv1yr, shape = rl)) + 
  scale_shape_manual(values = ss, breaks = c(1,0)) +
  scale_size(range = pp) +
  scale_color_manual(values = cc) +
  labs(y = NULL, x = NULL, title = spp[5], tag = '(e)') +
  mytheme + theme(plot.title = element_text(face = 'italic'))


## 2-YEAR SEEDLING SURVIVAL
### KP added ###
pred.sdlg$region <- as.character(pred.sdlg$region)
####
dat <- pred.sdlg %>% 
  mutate(prop18surv20 = if_else(is.na(prop18surv20), -1, prop18surv20)) %>% #give numerical values to NA values for point size
  mutate(region = if_else(prop18surv20 < 0 & region == 'MB', 'MB-NA', region)) %>% #different region name for NA for color
  mutate(region = if_else(prop18surv20 < 0 & region == 'RP', 'RP-NA', region)) %>% 
  mutate(rl = if_else(temprange_pos < -1, 0, 1)) %>% #create binary thermal range limit column
  mutate(rl = as.factor(rl)) #needs to be factor for shape 

plot6 <- dat %>%
  filter(full_species == spp[1]) %>% 
  ggplot(aes(x = elev, y = pred.2yr.surv, color = region)) + 
  geom_point(aes(size = prop18surv20, shape = rl)) + 
  scale_shape_manual(values = ss, breaks = c(1,0)) +
  scale_size(range = pp) +
  scale_color_manual(values = cc) +
  labs(y = NULL, x = NULL, title = spp[1], tag = '(f)') +
  mytheme + theme(plot.title = element_text(face = 'italic'))

plot7 <- dat %>%
  filter(full_species == spp[2]) %>% 
  ggplot(aes(x = elev, y = pred.2yr.surv, color = region)) + 
  geom_point(aes(size = prop18surv20, shape = rl)) + 
  scale_shape_manual(values = ss, breaks = c(1,0)) +
  scale_size(range = pp) +
  scale_color_manual(values = cc) +
  labs(y = NULL, x = NULL, title = spp[2], tag = '(g)') +
  mytheme + theme(plot.title = element_text(face = 'italic'))

plot8 <- dat %>%
  filter(full_species == spp[3]) %>% 
  ggplot(aes(x = elev, y = pred.2yr.surv, color = region)) + 
  geom_point(aes(size = prop18surv20, shape = rl)) + 
  scale_shape_manual(values = ss, breaks = c(1,0)) +
  scale_size(range = pp) +
  scale_color_manual(values = cc) +
  labs(y = NULL, x = NULL, title = spp[3], tag = '(h)') +
  mytheme + theme(plot.title = element_text(face = 'italic'))


plot9 <- dat %>%
  filter(full_species == spp[4]) %>% 
  ggplot(aes(x = elev, y = pred.2yr.surv, color = region)) + 
  geom_point(aes(size = prop18surv20, shape = rl)) + 
  scale_shape_manual(values = ss, breaks = c(1,0), labels = c("Below", "Beyond")) +
  scale_size(range = pp, guide = 'none') +
  scale_color_manual(values = cc, labels = c("West", "West-no recruitment", "East", "East-no recruitment")) +
  labs(y = NULL, x = NULL, title = spp[4], tag = '(i)') +
  mytheme + theme(plot.title = element_text(face = 'italic'), legend.position = 'bottom') +guides(color=guide_legend(title = 'Transect', title.position = 'top', nrow=1, byrow=TRUE,), shape = guide_legend(title = 'Range Limit Position', title.position = 'top', nrow = 1, byrow = TRUE)) 

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

plotLegend <- g_legend(plot9)
plot(plotLegend)

# Combine plots and save
fig <- grid.arrange(arrangeGrob(plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8, plot9 + 
                                  theme(legend.position = 'none'), heights=c(3,3,3)), 
                    plotLegend, heights = c(25,4), 
                    left = textGrob("Predicted Microhabitat Suitability\n(seedling survival)", #text for universal y-axis
                                    rot = 90, #rotate y-axis label 90 degrees
                                    gp = gpar(fontsize = tt)), #set font size for universal y-axis
                    bottom = textGrob('Elevation [m]', gp = gpar(fontsize = tt))
)

ggsave('outputs/ms_revisions_Dec2023/suitability-sdlgsurv.pdf', fig, width = 8, height = 8, units = c('in'))


####################################################################################################

# # FIG. S8 [MAT-ELEV] # # 

####################################################################################################

# Data
load('data/tidy/census_spp.RData') #adds relative recruitment to census.tms (microhab_analysis_prep.R)

dat <- census.spp %>% 
  distinct(census.spp$rep, .keep_all = TRUE) #keep only rep-level data

# Theme
pp <- 2 #point size
tt <- 12 #text size
cc <- c('#008000', 'black') #points color (https://html-color.codes/green)

mytheme <- theme(axis.title.x = element_text(size = tt), axis.text.x = element_text(colour = 'black', size = tt),
                 axis.title.y = element_text(size = tt), axis.text.y = element_text(colour = 'black', size = tt),
                 panel.background = element_blank(), axis.line = element_line(colour = "black"),
                 plot.tag = element_text(size = tt), plot.tag.position = c(0.035, 0.985),
                 legend.position = 'none')

## MAT ~ elev
# Test macroclim ~ elev relationship
mod <- lm(MAT ~ elev * region, data = dat) 
summary(mod) #only plot fitted lines if significant

# Plot results
MAT_elev <- ggplot(aes(x = elev, y = MAT, color = region), data = dat) +
  geom_point(size = pp) + scale_color_manual(values = cc) +
  labs(y = "Mean Annual Temperature [ºC]\n", x = NULL, tag = '(a)') +
  geom_smooth(method = "lm", se = FALSE) +
  mytheme


## MAP ~ elev
# Test macroclim ~ elev relationship
mod <- lm(MAP ~ elev * region, data = dat) 
summary(mod)

# Plot results
MAP_elev <- ggplot(aes(x = elev, y = MAP, color = region), data = dat) +
  geom_point(size = pp) + scale_color_manual(values = cc) +
  labs(y = "Mean Annual Precipitation [mm]\n", x = NULL, tag = '(b)') +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = cc, labels = c("West",  "East")) +
  mytheme + theme(legend.position = "none") + theme(plot.title = element_text(face = 'italic'), legend.position = 'bottom') + guides(color=guide_legend(title = 'Transect', title.position = 'top', nrow=1, byrow=TRUE,)) 

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

plotLegend <- g_legend(MAP_elev)

# Combine plots and save
fig <- grid.arrange(arrangeGrob(MAT_elev, MAP_elev + theme(legend.position = "none"), nrow = 1), plotLegend, heights = c(10,2), 
                    bottom = textGrob('Elevation [m]', gp = gpar(fontsize = tt)))

ggsave('outputs/ms_revisions_Dec2023/MAT_MAP_elev.pdf', fig, width = 8, height = 4, units = c('in'))

