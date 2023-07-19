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

unif_log_prior_try <- function(params, threshold, upper_mmax, b_value, epsilon, alpha = NULL, beta = NULL) {
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

gamma_unif_log_prior <- function(params, threshold, b_value, epsilon, alpha, beta, upper_mmax = NULL) {
  mmax <- params[1]
  mean <- params[2]
  stopifnot(!is.null(beta))
  stopifnot(!is.null(alpha))

  if((mmax <= mean) | (mean <= threshold) | (mean >= b_value+epsilon)){
    log_prior <- -1e7
  }
  else{
    log_prior <- -log(b_value + epsilon - threshold) + dgamma(mmax - mean, shape = alpha, rate = beta, log = TRUE)
  }

  return(log_prior)
}
