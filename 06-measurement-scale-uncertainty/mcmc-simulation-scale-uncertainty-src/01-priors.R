###########################
# LOG-PRIORS
###########################
# Import libraries
devtools::install_github('paulaoak/quaketools')
library(quaketools)

#Define the different uninformative prior functions for the GPD parameters


###########################
# CONSTANT THRESHOLD
###########################
unif_log_prior <- function(params) {
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

jeffreys_log_prior <- function(params) {
  sigma <- params[1]
  xi <- params[2]

  #jeffreys prior
  if((xi<-1/2) || (sigma<0)){
    log_prior <- -1e7
  }
  else{
    log_prior <- -log(sigma)-0.5*log(2*xi+1)-log(xi+1)
  }

  return(log_prior)
}

mdi_log_prior <- function(params) {
  sigma <- params[1]
  xi <- params[2]

  #mdi prior
  if((xi<-1) || (sigma<0)){
    log_prior <- -1e7
  }
  else{
    log_prior <- -log(sigma)-xi
  }

  return(log_prior)
}

