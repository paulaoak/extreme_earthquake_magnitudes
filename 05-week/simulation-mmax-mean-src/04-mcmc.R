#################################
# MH WITHIN GIBBS ALGORITHM
#################################

# Perform MCMC sampling using the Metropolis-Hastings within Gibbs algorithm

###########################
# CONSTANT THRESHOLD UNROUNDED OBSERVATIONS
###########################
bayesian_estimation_gpd_mmax_mean <- function(mmax_mean_initial, x, u, n_iter, n_burn, sd_mmax, sd_mean, prior_choice = c('unif-unif', 'unif-gamma', 'flat-flat', 'flat-gamma', 'gamma-gamma'), b_value = NULL, epsilon = NULL, alpha1 = NULL, beta1 = NULL, upper_mmax = NULL, alpha2 = NULL, beta2 = NULL){

  prior_choice <- match.arg(prior_choice)
  if(prior_choice == 'unif-unif'){
    prior_choice <- unif_log_prior
  }
  else if(prior_choice == 'unif-gamma'){
    prior_choice <- gamma_unif_log_prior
  }
  else if(prior_choice == 'flat-flat'){
    prior_choice <- flat_log_prior
  }
  else if(prior_choice == 'flat-gamma'){
    prior_choice <- gamma_flat_log_prior
  }
  else if(prior_choice == 'gamma-gamma'){
    prior_choice <- gamma_gamma_log_prior
  }

  # Initial value of parameters
  mmax <- mmax_mean_initial[1]
  mean <- mmax_mean_initial[2]

  # Storage for sampled parameter values
  params_samples <- matrix(NA, ncol = n_iter, nrow = 2)

  # Metropolis-Hastings within Gibbs algorithm
  for (iter in 1:n_iter) {

    #MH for mmax
    min_proposal <- u
    log_posterior_mmax <- function (mmax_1){log_posterior_gpd_mmax_mean(params = c(mmax_1, mean), x = x, threshold = u, prior = prior_choice, upper_mmax = upper_mmax, b_value = b_value, epsilon = epsilon, alpha1 = alpha1, beta1 = beta1, alpha2 = alpha2, beta2 = beta2)}
    mmax <- mh_step_random_walk_positive (mmax, log_posterior_mmax, sd = sd_mmax, a = min_proposal, b = 15)

    #MH for mean
    log_posterior_mean <- function (mean_1){log_posterior_gpd_mmax_mean(params = c(mmax, mean_1), x = x, threshold = u, prior = prior_choice, upper_mmax = upper_mmax, b_value = b_value, epsilon = epsilon, alpha = alpha, beta = beta, alpha2 = alpha2, beta2 = beta2)}
    mean <- mh_step_random_walk_positive (mean, log_posterior_mean, sd = sd_mean, a = u, b = 15)
    print(c(mmax, mean, log_posterior_mmax(mmax), log_posterior_mean(mean)))

    params_samples[,iter] <- c(mmax, mean)
  }

  params_samples[,n_burn:n_iter]
}

