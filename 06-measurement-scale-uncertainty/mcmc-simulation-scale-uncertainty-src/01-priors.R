###########################
# LOG-PRIORS
###########################
# Import libraries
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

problem_context_log_prior <- function(params) {
  sigma <- params[1]
  xi <- params[2]

  #problem context prior
  if((xi<-1) || (sigma<0) || (xi>1)){
    log_prior <- -1e7
  }
  else{
    xi_mod <- (1+xi)/2
    logit_xi <- log(xi_mod/(1-xi_mod))
    prior <- dnorm(logit_xi, mean = 0.42, sd = 0.4) * dgamma(sigma, 9, 10)
    log_prior <- log(prior)
  }

  return(log_prior)
}
