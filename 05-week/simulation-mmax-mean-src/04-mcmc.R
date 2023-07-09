#################################
# MH WITHIN GIBBS ALGORITHM
#################################

# Perform MCMC sampling using the Metropolis-Hastings within Gibbs algorithm

###########################
# CONSTANT THRESHOLD UNROUNDED OBSERVATIONS
###########################
bayesian_estimation_gpd_mmax_mean <- function(mmax_mean_initial, x, u, n_iter, n_burn, sd_mmax, sd_mean, prior_choice = c('flat-flat', 'flat-gamma'), b_value, epsilon, alpha = NULL, beta = NULL, upper_mmax = NULL){

  prior_choice <- match.arg(prior_choice)
  if(prior_choice == 'flat-flat'){
    prior_choice <- unif_log_prior
  }
  else if(prior_choice == 'flat-gamma'){
    prior_choice <- unif_log_prior
  }
  else{
    prior_choice <- gamma_unif_log_prior
  }

  # Initial value of parameters
  mmax <- mmax_mean_initial[1]
  mean <- mmax_mean_initial[2]

  # Storage for sampled parameter values
  params_samples <- matrix(NA, ncol = n_iter, nrow = 2)

  # Metropolis-Hastings within Gibbs algorithm
  for (iter in 1:n_iter) {

    #MH for scale parameter
    log_posterior_mmax <- function (mmax_1){log_posterior_gpd_mmax_mean(params = c(mmax_1, mean), x = x, threshold = u, prior = prior_choice, upper_mmax = upper_mmax, b_value = b_value, epsilon = epsilon, alpha = alpha, beta = beta)}
    mmax <- mh_step_random_walk_positive (mmax, log_posterior_mmax, sd = sd_mmax)

    #MH for shape parameter
    log_posterior_mean <- function (mean_1){log_posterior_gpd_mmax_mean(params = c(mmax, mean_1), x = x, threshold = u, prior = prior_choice, upper_mmax = upper_mmax, b_value = b_value, epsilon = epsilon, alpha = alpha, beta = beta)}
    mean <- mh_step_random_walk_positive (mean, log_posterior_mean, sd = sd_mean)

    params_samples[,iter] <- c(mmax, mean)
  }

  params_samples[,n_burn:n_iter]
}

