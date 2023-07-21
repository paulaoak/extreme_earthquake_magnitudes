###########################
# LOG-PRIORS
###########################

# Import libraries
library(quaketools)

#Define the different uninformative prior functions for the Mmax and mean
#parameters of a GPD distribution


###########################
# CONSTANT THRESHOLD
###########################

#Uniform prior on mean parameter between threshold and b_value approximation +
#an epsilon error and the mmax parameter uniform prior between the mean and an
#upper bound for mmax
#It is important to notice that in this case mmax and mean are not independent
#unlike when we specified a prior on the shape and scale parameters

unif_log_prior <- function(params, threshold, upper_mmax, b_value, epsilon, alpha1 = NULL, beta1 = NULL, alpha2=NULL, beta2=NULL) {
  mmax <- params[1]
  mean <- params[2]
  stopifnot(!is.null(upper_mmax))

  stopifnot(b_value + epsilon < upper_mmax)

  if((mean >= mmax) | (mean <= threshold) | (mean >= b_value+epsilon) | (mmax >= upper_mmax)){
    log_prior <- -1e7
  }
  else{
    prior <- 1/((upper_mmax - threshold)^2 - (upper_mmax - b_value - epsilon)^2)
    log_prior <- log(prior)
  }

  return(log_prior)
}

#Uniform prior on mean parameter between threshold and b_value approximation +
#an epsilon error and the difference mmax-mean has a gamma distribution
#As before, mmax and mean are not independent

gamma_unif_log_prior <- function(params, threshold, b_value, epsilon, alpha1, beta1, upper_mmax = NULL, alpha2=NULL, beta2=NULL) {
  mmax <- params[1]
  mean <- params[2]
  stopifnot(!is.null(beta1))
  stopifnot(!is.null(alpha1))

  if((mmax <= mean) | (mean <= threshold) | (mean >= b_value+epsilon)){
    log_prior <- -1e7
  }
  else{
    log_prior <- -log(b_value + epsilon - threshold) + dgamma(mmax - mean, shape = alpha1, rate = beta1, log = TRUE)
  }

  return(log_prior)
}

#Flat uninformative prior on mean and Mmax only including the constraint
#threshold<mean<mmax
flat_log_prior <- function(params, threshold, upper_mmax = NULL, b_value= NULL, epsilon= NULL, alpha1 = NULL, beta1 = NULL, alpha2=NULL, beta2=NULL) {
  mmax <- params[1]
  mean <- params[2]

  if((mmax <= mean) | (mean <= threshold)){
    log_prior <- -1e7
  }
  else{
    log_prior <- 0
  }

  return(log_prior)
}

#Improper flat prior on the mean and gamma prior on the difference Mmax-mean
gamma_flat_log_prior <- function(params, threshold, b_value = NULL, epsilon = NULL, alpha1, beta1, upper_mmax = NULL, alpha2=NULL, beta2=NULL) {
  mmax <- params[1]
  mean <- params[2]
  stopifnot(!is.null(beta1))
  stopifnot(!is.null(alpha1))

  if((mmax <= mean) | (mean <= threshold)){
    log_prior <- -1e7
  }
  else{
    log_prior <- dgamma(mmax - mean, shape = alpha1, rate = beta1, log = TRUE)
  }

  return(log_prior)
}


#Improper flat prior on the mean and gamma prior on the difference Mmax-mean
gamma_gamma_log_prior <- function(params, threshold, alpha1, beta1, alpha2, beta2, upper_mmax=NULL, b_value=NULL, epsilon=NULL) {
  mmax <- params[1]
  mean <- params[2]
  stopifnot(!is.null(beta1))
  stopifnot(!is.null(alpha1))
  stopifnot(!is.null(beta2))
  stopifnot(!is.null(alpha2))

  if((mmax <= mean) | (mean <= threshold)){
    log_prior <- -1e7
  }
  else{
    log_prior <- dgamma(mean - threshold, shape = alpha1, rate = beta1, log = TRUE) +
                 dgamma(mmax - mean, shape = alpha2, rate = beta2, log = TRUE)
  }

  return(log_prior)
}

