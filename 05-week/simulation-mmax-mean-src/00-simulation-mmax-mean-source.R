# Import necessary libraries
library(quaketools)

# Source necessary functions for the analysis
source(here::here("05-week", "simulation-mmax-mean-src",
                  "01-priors.R"))

source(here::here("05-week", "simulation-mmax-mean-src",
                  "02-posteriors.R"))

source(here::here("05-week", "simulation-mmax-mean-src",
                  "03-mh-step.R"))

source(here::here("05-week", "simulation-mmax-mean-src",
                  "04-mcmc.R"))
