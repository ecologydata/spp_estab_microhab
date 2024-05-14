last updated: 14 May 2024

# Microhabitat Analyses

Description of scripts and data for Chardon et al. 2024. Variable species establishment in response to microhabitat indicates different likelihoods of climate-driven range shifts. Ecography. Files are ordered under their parent folder.


## Compiled data

data/abbreviations.xlsx: description of all variables (not all used in scripts below) and species abbreviations

data/census_sdlg.RData: yearly seedling survival + microhabitat variables

data/census_spp.RData: yearly census data + microhabitat variables at quadrat level

data/census_rep.RData: yearly census data + microhabitat variables at rep level

data/thermal2elevRL.RData: thermal to elev range limit per species

data/seedmix_added: weight or number of seeds added to seedmix per species


## Analyses (results saved to 'output/')

scripts/01_data_expl.R: Choose best random effects structure, variance partitioning, choose fixed effects

scripts/02_param_corr.R: choose biologically meaningful explanatory parameters for models

scripts/03_varsel_recruitprob_bin.R: Binomial GLM for recruitment probability

scripts/04_varsel_recruit_bin.R: Binomial GLM for relative number of recruits

scripts/05_varsel_sdlg_bin.R: Binomial GLM for 1- and 2-year seedling survival

scripts/06_therm2elevRL.R: calculate elevational range limit corresponding to thermal RL

scripts/07_microhab_suitab.R: predict binary/continuous recruitment and 1/2 year survival to each plot


## Results (tables saved to 'output/' and figure files not saved)

scripts/ms_figs.R: figures for manuscript (note: 'figs/' folder is set up to run this script and save figures)

scripts/ms_tables.R: tables for manuscript