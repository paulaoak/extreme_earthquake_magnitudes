###########################
# CONSTANT THRESHOLD UNROUNDED OBSERVATIONS
###########################
bayesian_estimation_gpd <- function(sigxi_initial, x, u, n_iter, n_burn, sd_xi, sd_sig_u, prior_choice = c('flat', 'mdi', 'jeffreys')){

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
  sig_u <- sigxi_initial[1]
  xi <- sigxi_initial[2]

  # Storage for sampled parameter values
  params_samples <- matrix(NA, ncol = n_iter, nrow = 2)

  # Metropolis-Hastings within Gibbs algorithm
  for (iter in 1:n_iter) {

    #MH for scale parameter
    log_posterior_sig_u <- function (scale){log_posterior_gpd(params = c(scale, xi), x = x, u = u, prior = prior_choice)}
    sig_u <- mh_step_random_walk_positive (sig_u, log_posterior_sig_u, sd = sd_sig_u, a = 0)

    #MH for shape parameter
    log_posterior_xi <- function (shape){log_posterior_gpd(params = c(sig_u, shape), x = x, u = u, prior = prior_choice)}
    xi <- mh_step_random_walk (xi, log_posterior_xi, sd = sd_xi)

    params_samples[,iter] <- c(sig_u, xi)
  }

  params_samples[,n_burn:n_iter]
}
