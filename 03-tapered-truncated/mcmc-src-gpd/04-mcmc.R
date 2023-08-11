###########################
# CONSTANT THRESHOLD UNROUNDED OBSERVATIONS NO SCALE
###########################
bayesian_estimation_gpd <- function(sigxi_initial, x, u, n_iter, n_burn, sd_xi, sd_sig_u, prior_choice = c('flat')){

  prior_choice <- match.arg(prior_choice)
  if(prior_choice == 'flat'){
    prior_choice <- flat_log_prior
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


# NOT TAKING INTO ACCOUNT SCALE UNCERTAINTY
simulation_mcmc_gpd <- function(xi, sigma,
                            xi_init_min = -0.5, sigma_init_min = 0,
                            xi_init_max = 1.5, sigma_init_max = 4,
                            n_starting_points = 5,
                            u,
                            n_data = 500, n_iter = 1e4, n_burn = 1e3,
                            prior = c('flat'),
                            sd_xi = c(0.2, 0.15, 0.1, 0.05), sd_sig = c(0.2, 0.15, 0.1, 0.05)){

  prior <- match.arg(prior)

  # Generate data
  set.seed(123)
  x <- rgpd(n = n_data, scale = sigma, shape = xi, shift = u)

  #Perform grid search to find optimal variances for the random walk
  n_iter_grid <- 1e3
  n_burn_grid <- 1e2
  grid <- crossing(sd_xi, sd_sig)
  grid <- as.data.frame(grid)
  n_combinations <- dim(grid)[1]
  set.seed(111)
  params_initial_grid <- matrix(c(runif(n_starting_points, min = sigma_init_min, max = sigma_init_max),
                                  runif(n_starting_points, min = xi_init_min, max = xi_init_max)),
                                ncol = 2)

  ess_grid <- sapply(1:n_combinations, function(i){
    sd_xi <- grid[i, 1]
    sd_sig <- grid[i, 2]
    ess_start <- sapply(1:n_starting_points, function(j){

      params_samples <- bayesian_estimation_gpd(sigxi_initial = params_initial_grid[j,], x = x, u = u, n_iter = n_iter_grid,
                                                n_burn = n_burn_grid, sd_xi = sd_xi, sd_sig = sd_sig, prior_choice = prior)
      effectiveSize(t(params_samples))

    })
    mean(ess_start)
  })
  optim_sig_xi_sig <- grid[which.max(ess_grid),]

  # Run mcmc chains with the optimal variances obtained for different starting
  # points and asses convergence
  posterior_samples <- matrix(NA, nrow = 2, ncol = (n_iter - n_burn + 1) * n_starting_points) #preallocate space
  posterior_sigma <- list()
  posterior_xi <- list()



  for(j in 1:n_starting_points) {
    posterior_samples[1:2, ((n_iter - n_burn + 1) * (j - 1) + 1): ((n_iter - n_burn + 1) * j)] <- bayesian_estimation_gpd(sigxi_initial = params_initial_grid[j,], x = x, u = u, n_iter = n_iter,
                                                                                                                          n_burn = n_burn, sd_xi = optim_sig_xi_sig[[1]], sd_sig = optim_sig_xi_sig[[2]], prior_choice = prior)
    posterior_sigma[[j]] <- mcmc(posterior_samples[1, ((n_iter - n_burn + 1) * (j - 1) + 1): ((n_iter - n_burn + 1) * j)])
    posterior_xi[[j]] <- mcmc(posterior_samples[2, ((n_iter - n_burn + 1) * (j - 1) + 1): ((n_iter - n_burn + 1) * j)])
  }

  # Convergence diagnostics
  upper_ci_sigma <- gelman.diag(mcmc.list(posterior_sigma))[1]$psrf[2]
  #stopifnot('Chains with samples of sigma have not converge to the same distribution.' = upper_ci_sigma < 1.2)
  upper_ci_xi <- gelman.diag(mcmc.list(posterior_xi))[1]$psrf[2]
  #stopifnot('Chains with samples of xi have not converge to the same distribution.' = upper_ci_xi < 1.2)

  #if((upper_ci_sigma_u > 1.25) || (upper_ci_xi > 1.25)){posterior_samples<- c(0,0)}
  print(c(upper_ci_sigma, upper_ci_xi))
  return(posterior_samples)

}

