# MH WITHIN GIBBS ALGORITHM-----------------------------------------------------
# Perform MCMC sampling using the Metropolis-Hastings within Gibbs algorithm
bayesian_estimation_tgpd <- function(betatheta_initial, x, u, n_iter, n_burn, sd_beta, sd_theta){
  # Initial value of parameters
  beta <- betatheta_initial[1]
  theta <- betatheta_initial[2]

  # Storage for sampled parameter values
  params_samples <- matrix(NA, ncol = n_iter, nrow = 2)

  # Metropolis-Hastings within Gibbs algorithm
  for (iter in 1:n_iter) {

    #MH for scale parameter
    log_posterior_beta <- function (power){log_posterior_tapered(params = c(power, theta), x = x, u = u, prior = unif_log_prior)}
    beta <- mh_step_random_walk_positive (beta, log_posterior_beta, sd = sd_beta, a = 0)

    #MH for taper scale parameter
    log_posterior_theta <- function (scale_taper){log_posterior_tapered(params = c(beta, scale_taper), x = x, u = u, prior = unif_log_prior)}
    theta <- mh_step_random_walk_positive (theta, log_posterior_theta, sd = sd_theta, a = 0)

    params_samples[,iter] <- c(beta, theta)
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

simulation_mcmc_tapered <- function(xi, sigma,
                                    beta_init_min = 0, theta_init_min = 0,
                                    beta_init_max = 3, theta_init_max = 4,
                                    n_starting_points = 5, u,
                                    n_data, n_iter = 1e4, n_burn = 1e3,
                                    sd_beta = c(0.2, 0.15, 0.1, 0.05),
                                    sd_theta = c(0.5, 0.2, 0.15, 0.1, 0.05)){


  # Generate GPD data
  set.seed(123)
  x <- rgpd(n = n_data, scale = sigma, shape = xi, shift = u)

  #Perform grid search to find optimal variances for the random walk
  n_iter_grid <- 1e3
  n_burn_grid <- 1e2
  grid <- crossing(sd_beta, sd_theta)
  grid <- as.data.frame(grid)
  n_combinations <- dim(grid)[1]
  set.seed(111)
  params_initial_grid <- matrix(c(runif(n_starting_points, min = beta_init_min, max = beta_init_max),
                                  runif(n_starting_points, min = theta_init_min, max = theta_init_max)),
                                  ncol = 2)

  ess_grid <- sapply(1:n_combinations, function(i){
    sd_beta <- grid[i, 1]
    sd_theta <- grid[i, 2]
    ess_start <- sapply(1:n_starting_points, function(j){
      params_samples <- bayesian_estimation_tgpd(betatheta_initial = params_initial_grid[j,], x = x, u = u, n_iter = n_iter_grid,
                                                 n_burn = n_burn_grid, sd_beta = sd_beta, sd_theta = sd_theta)
      effectiveSize(t(params_samples))

    })
    mean(ess_start)
  })
  optim_sig_xi_sig <- grid[which.max(ess_grid),]

  # Run mcmc chains with the optimal variances obtained for different starting
  # points and asses convergence
  posterior_samples <- matrix(NA, nrow = 2, ncol = (n_iter - n_burn + 1) * n_starting_points) #preallocate space
  posterior_beta <- list()
  posterior_theta <- list()

  for(j in 1:n_starting_points) {
    posterior_samples[, ((n_iter - n_burn + 1) * (j - 1) + 1): ((n_iter - n_burn + 1) * j)] <- bayesian_estimation_tgpd(betatheta_initial = params_initial_grid[j,], x = x, u = u, n_iter = n_iter,
                                                                                                                        n_burn = n_burn, sd_beta = optim_sig_xi_sig[[1]], sd_theta = optim_sig_xi_sig[[2]])

    posterior_beta[[j]] <- mcmc(posterior_samples[1, ((n_iter - n_burn + 1) * (j - 1) + 1): ((n_iter - n_burn + 1) * j)])
    posterior_theta[[j]] <- mcmc(posterior_samples[2, ((n_iter - n_burn + 1) * (j - 1) + 1): ((n_iter - n_burn + 1) * j)])
  }

  # Convergence diagnostics
  upper_ci_beta <- gelman.diag(mcmc.list(posterior_beta))[1]$psrf[2]
  #stopifnot('Chains with samples of sigma_u have not converge to the same distribution.' = upper_ci_sigma_u < 1.2)
  upper_ci_theta <- gelman.diag(mcmc.list(posterior_theta))[1]$psrf[2]

  #if((upper_ci_sigma > 1.25) || (upper_ci_xi > 1.25) || (upper_ci_lambda > 1.25)){posterior_samples<- c(0,0)}
  print(c(upper_ci_beta, upper_ci_theta))

  return(posterior_samples)
}
