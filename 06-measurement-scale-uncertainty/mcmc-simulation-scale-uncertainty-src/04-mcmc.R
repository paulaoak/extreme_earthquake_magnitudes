#################################
# MH WITHIN GIBBS ALGORITHM
#################################

# Perform MCMC sampling using the Metropolis-Hastings within Gibbs algorithm

###########################
# GPD ACCOUNTING FOR ROUNDING UNDER PENULTIMATE APPROXIMATION
###########################
bayesian_estimation_gpd_scale_penultimate <- function(sigxi_lambda_initial, c, x, u, n_iter, n_burn, sd_xi, sd_sigma, min_lambda, max_lambda, prior_choice = c('flat', 'mdi', 'jeffreys')){

  prior_choice <- match.arg(prior_choice)
  if(prior_choice == 'flat'){
    prior_choice <- unif_log_prior
  }
  else if(prior_choice == 'mdi'){
    prior_choice <- mdi_log_prior
  }
  else{
    prior_choice <- jeffreys_log_prior
  }

  # Initial value of parameters
  sig <- sigxi_lambda_initial[1]
  xi <- sigxi_lambda_initial[2]
  lambda <- sigxi_lambda_initial[3]

  # Storage for sampled parameter values
  params_samples <- matrix(NA, ncol = n_iter, nrow = 3)

  # Metropolis-Hastings within Gibbs algorithm
  for (iter in 1:n_iter) {

    #MH for scale parameter
    log_posterior_sig <- function (scale){log_posterior_gpd_scale_penultimate(params = c(scale, xi, lambda), c = c, x = x, u = u, prior = prior_choice)}
    sig <- mh_step_random_walk_positive (sig, log_posterior_sig, sd = sd_sigma, a = 0)

    #MH for shape parameter
    log_posterior_xi <- function (shape){log_posterior_gpd_scale_penultimate(params = c(sig, shape, lambda), c = c, x = x, u = u, prior = prior_choice)}
    xi <- mh_step_random_walk (xi, log_posterior_xi, sd = sd_xi)

    #MH for lambda parameter
    lambda <- mh_step_independent (min_proposal = min_lambda, max_proposal = max_lambda)

    params_samples[,iter] <- c(sig, xi, lambda)

  }

  params_samples[,n_burn:n_iter]
}


###########################
# GPD ACCOUNTING FOR ROUNDING UNDER ASSUMPTION OF VALIDITY OF GPD FOR DIFFERENT LAMBDAS
###########################
bayesian_estimation_gpd_scale <- function(sigxi_lambda_initial, x, u, n_iter, n_burn, sd_xi, sd_sigma, min_lambda, max_lambda, prior_choice = c('flat', 'mdi', 'jeffreys')){

  prior_choice <- match.arg(prior_choice)
  if(prior_choice == 'flat'){
    prior_choice <- unif_log_prior
  }
  else if(prior_choice == 'mdi'){
    prior_choice <- mdi_log_prior
  }
  else{
    prior_choice <- jeffreys_log_prior
  }

  # Initial value of parameters
  sig <- sigxi_lambda_initial[1]
  xi <- sigxi_lambda_initial[2]
  lambda <- sigxi_lambda_initial[3]

  # Storage for sampled parameter values
  params_samples <- matrix(NA, ncol = n_iter, nrow = 3)

  # Metropolis-Hastings within Gibbs algorithm
  for (iter in 1:n_iter) {

    #MH for scale parameter
    log_posterior_sig <- function (scale){log_posterior_gpd_scale(params = c(scale, xi, lambda), x = x, u = u, prior = prior_choice)}
    sig <- mh_step_random_walk_positive (sig, log_posterior_sig, sd = sd_sigma, a = 0)

    #MH for shape parameter
    log_posterior_xi <- function (shape){log_posterior_gpd_scale(params = c(sig, shape, lambda), x = x, u = u, prior = prior_choice)}
    xi <- mh_step_random_walk (xi, log_posterior_xi, sd = sd_xi)

    #MH for lambda parameter
    lambda <- mh_step_independent (min_proposal = min_lambda, max_proposal = max_lambda)

    params_samples[,iter] <- c(sig, xi, lambda)

  }

  params_samples[,n_burn:n_iter]
}
