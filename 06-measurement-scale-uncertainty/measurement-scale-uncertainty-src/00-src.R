# Import necessary libraries
library(quaketools)
library(cooltools)

# Source necessary functions for the analysis
source(here::here("06-measurement-scale-uncertainty", "measurement-scale-uncertainty-src",
                  "01-likelihood-scale-uncertainty.R"))

source(here::here("06-measurement-scale-uncertainty", "measurement-scale-uncertainty-src",
                  "02-mle-scale-uncertainty.R"))

source(here::here("06-measurement-scale-uncertainty", "measurement-scale-uncertainty-src",
                  "03-profile-likelihood-scale.R"))

#source(here::here("06-measurement-scale-uncertainty", "measurement-scale-uncertainty-src",
#                  "04-contours-profile-likelihood.R"))

source(here::here("06-measurement-scale-uncertainty", "measurement-scale-uncertainty-src",
                  "05-weighted-least-squares.R"))
