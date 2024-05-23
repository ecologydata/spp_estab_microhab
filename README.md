last updated: 23 May 2024

## Description of the data and file structure

Made available here are compiled data on yearly germination/seedling census and microhabitat variables, and R scripts used to analyze the effects of microhabitat on germination and seedling survival. These data and R scripts were used in Chardon et al. 2024. Variable species establishment in response to microhabitat indicates different likelihoods of climate-driven range shifts. Ecography. Files are listed under their parent folder.

## Compiled Data: these data files are used in the R scripts below, which load these data files as listed

data/abbreviations.xlsx: description of all variables and species abbreviations, with datasets each variable is (1) or is not (0) used in

data/census_sdlg.RData: yearly seedling survival + microhabitat variables

data/census_spp.RData: yearly census data + microhabitat variables at quadrat level (3 quadrats per rep)

data/census_rep.RData: yearly census data + microhabitat variables at rep level (3 reps per block, with 2 block per site)

data/thermal2elevRL.RData: thermal to elevational range limit per species

data/seedmix_added: weight or number of seeds added to seedmix per species

## Analyses: to be run in numeric order (uploaded to Zenodo)

scripts/01_data_expl.R: choose best random effects structure, variance partitioning, choose fixed effects

scripts/02_param_corr.R: choose biologically meaningful explanatory parameters for models

scripts/03_varsel_recruitprob_bin.R: Binomial GLM for recruitment probability

scripts/04_varsel_recruit_bin.R: Binomial GLM for relative number of recruits

scripts/05_varsel_sdlg_bin.R: Binomial GLM for 1- and 2-year seedling survival

scripts/06_therm2elevRL.R: calculate elevational range limit corresponding to thermal RL

scripts/07_microhab_suitab.R: predict binary/continuous recruitment and 1/2 year survival to each plot

## Results: scripts used to generate figures and tables for manuscript (uploaded to Zenodo)

scripts/ms_figs.R: figures for manuscript

scripts/ms_tables.R: tables for manuscript

## Package versioning

The folder renv/ contains the packages versions used to run the above scripts. To use, download the entire repository and run renv::update().

## Sharing/Access information

Available as GitHub repository: https://github.com/nchardon/spp_estab_microhab/

## Code/Software

We conducted all data processing and analyses in R version 4.2.3 (R Core Team 2023).

## References

R Core Team, R. 2023. “R: A Language and Environment for Statistical Computing.” R Foundation for Statistical Computing, Vienna, Austria. http://www.r-project.org.