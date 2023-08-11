###########################
# LOG-PRIORS
###########################
# Import libraries
library(quaketools)

#Define the different uninformative prior functions for the GPD parameters


###########################
# CONSTANT THRESHOLD
###########################
flat_log_prior <- function(params) {
  sigma <- params[1]

  #uniform prior for xi and log(sigma)
  if(sigma<0){
    log_prior <- -1e7
  }
  else{
    log_prior <- -log(sigma)
  }

  return(log_prior)
}
