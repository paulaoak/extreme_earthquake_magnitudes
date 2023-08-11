# MH WITHIN GIBBS ALGORITHM-----------------------------------------------------
# Perform MCMC sampling using the Metropolis-Hastings within Gibbs algorithm
bayesian_estimation_truncated_gr <- function(beta_mmax_initial, x, u, n_iter, n_burn, sd_beta, sd_mmax){
  # Initial value of parameters
  beta <- beta_mmax_initial[1]
  mmax <- beta_mmax_initial[2]

  # Storage for sampled parameter values
  params_samples <- matrix(NA, ncol = n_iter, nrow = 2)

  # Metropolis-Hastings within Gibbs algorithm
  for (iter in 1:n_iter) {

    #MH for scale parameter
    log_posterior_beta <- function (rate){log_posterior_truncated_gr(params = c(rate, mmax), x = x, u = u, prior = unif_log_prior_gr)}
    beta <- mh_step_random_walk_positive (beta, log_posterior_beta, sd = sd_beta, a = 0)

    #MH for shape parameter
    log_posterior_mmax <- function (maximum){log_posterior_truncated_gr(params = c(beta, maximum), u = u, x = x, prior = unif_log_prior_gr)}
    mmax <- mh_step_random_walk_positive (mmax, log_posterior_mmax, sd = sd_mmax, a = u)

    params_samples[,iter] <- c(beta, mmax)
  }

  params_samples[,n_burn:n_iter]
}

#################################
##Function for simulation
#################################
###############
#SIMULATIONS
###############

# Import libraries
library(coda)
library(tidyr)

simulation_mcmc_truncated_gr <- function(xi, sigma,
                                         beta_init_min = 0,
                                         beta_init_max = 4, mmax_init_max = 10,
                                         n_starting_points = 5, u,
                                         n_data, n_iter = 1e4, n_burn = 1e3,
                                         sd_beta = c(0.2, 0.15, 0.1),
                                         sd_mmax = c(2.25, 1.50, 0.75, 0.2)){


  # Generate GPD data
  set.seed(123)
  x <- rgpd(n = n_data, scale = sigma, shape = xi, shift = u)

  #Perform grid search to find optimal variances for the random walk
  n_iter_grid <- 1e3
  n_burn_grid <- 1e2
  grid <- crossing(sd_beta, sd_mmax)
  grid <- as.data.frame(grid)
  n_combinations <- dim(grid)[1]
  set.seed(111)
  params_initial_grid <- matrix(c(runif(n_starting_points, min = beta_init_min, max = beta_init_max),
                                  runif(n_starting_points, min = u, max = mmax_init_max)),
                                  ncol = 2)

  ess_grid <- sapply(1:n_combinations, function(i){
    sd_beta <- grid[i, 1]
    sd_mmax <- grid[i, 2]

    ess_start <- sapply(1:n_starting_points, function(j){
      params_samples <- bayesian_estimation_truncated_gr(beta_mmax_initial = params_initial_grid[j,], x = x, u = u, n_iter = n_iter_grid,
                                                         n_burn = n_burn_grid, sd_beta = sd_beta, sd_mmax = sd_mmax)
      effectiveSize(t(params_samples))

    })
    mean(ess_start)
  })
  optim_sig_xi_sig <- grid[which.max(ess_grid),]

  # Run mcmc chains with the optimal variances obtained for different starting
  # points and asses convergence
  posterior_samples <- matrix(NA, nrow = 2, ncol = (n_iter - n_burn + 1) * n_starting_points) #preallocate space
  posterior_beta <- list()
  posterior_mmax <- list()


  for(j in 1:n_starting_points) {
    posterior_samples[, ((n_iter - n_burn + 1) * (j - 1) + 1): ((n_iter - n_burn + 1) * j)] <- bayesian_estimation_truncated_gr(beta_mmax_initial = params_initial_grid[j,], x = x, u = u, n_iter = n_iter,
                                                                                                                                n_burn = n_burn, sd_beta = optim_sig_xi_sig[[1]], sd_mmax = optim_sig_xi_sig[[2]])

    posterior_beta[[j]] <- mcmc(posterior_samples[1, ((n_iter - n_burn + 1) * (j - 1) + 1): ((n_iter - n_burn + 1) * j)])
    posterior_mmax[[j]] <- mcmc(posterior_samples[2, ((n_iter - n_burn + 1) * (j - 1) + 1): ((n_iter - n_burn + 1) * j)])

  }

  # Convergence diagnostics
  upper_ci_beta <- gelman.diag(mcmc.list(posterior_beta))[1]$psrf[2]
  #stopifnot('Chains with samples of sigma_u have not converge to the same distribution.' = upper_ci_sigma_u < 1.2)
  upper_ci_mmax <- gelman.diag(mcmc.list(posterior_mmax))[1]$psrf[2]
  #stopifnot('Chains with samples of xi have not converge to the same distribution.' = upper_ci_xi < 1.2)

  #if((upper_ci_sigma > 1.25) || (upper_ci_xi > 1.25) || (upper_ci_lambda > 1.25)){posterior_samples<- c(0,0)}
  print(c(upper_ci_beta, upper_ci_mmax))

  return(posterior_samples)
}
