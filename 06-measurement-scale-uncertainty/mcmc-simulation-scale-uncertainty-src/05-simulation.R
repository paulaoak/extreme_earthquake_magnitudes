###############
#SIMULATIONS
###############

## PENULTIMATE APPROXIMATION

# Import libraries
library(coda)
library(tidyr)

simulation_mcmc_scale_penultimate <- function(xi, sigma,
                                              xi_init_min = -0.5, sigma_u_init_min = 0,
                                              xi_init_max = 1.5, sigma_u_init_max = 4,
                                              n_starting_points = 5, u,
                                              n_data = 500, n_iter = 1e4, n_burn = 1e3,
                                              prior = c('mdi', 'flat', 'jeffreys'),
                                              sd_xi = c(0.2, 0.15, 0.1, 0.05), sd_sigma = c(0.2, 0.15, 0.1, 0.05),
                                              min_lambda = -0.1, max_lambda = 3, step_lambda = 0.1,
                                              min_xi = -0.35, max_xi = 0.35, step_xi = 0.05,
                                              scale_transformation = TRUE){

  prior <- match.arg(prior)

  # Obtain estimate of c using weighted least squares regression
  wlm_c <- estimate_c_penultimate(min_lambda = min_lambda, max_lambda = max_lambda, step_lambda = step_lambda,
                                  min_xi = min_xi, max_xi = max_xi, step_xi = step_xi, u = u, x = x)
  c <- wlm_c[1]

  # Generate data
  set.seed(123)
  x <- rgpd(n = n_data, scale = sigma, shape = xi, shift = u)
  if (scale_transformation == TRUE){
    # Data under scale transformation
    x <- log(x)
  }

  #Perform grid search to find optimal variances for the random walk
  n_iter_grid <- 1e3
  n_burn_grid <- 1e2
  grid <- crossing(sd_xi, sd_sigma)
  grid <- as.data.frame(grid)
  n_combinations <- dim(grid)[1]
  set.seed(111)
  params_initial_grid <- matrix(c(runif(n_starting_points, min = sigma_u_init_min, max = sigma_u_init_max),
                                  runif(n_starting_points, min = xi_init_min, max = xi_init_max),
                                  runif(n_starting_points, min = min_lambda, max = max_lambda)),
                                ncol = 3)

  ess_grid <- sapply(1:n_combinations, function(i){
    sd_xi <- grid[i, 1]
    sd_sigma <- grid[i, 2]
    ess_start <- sapply(1:n_starting_points, function(j){

      params_samples <- bayesian_estimation_gpd_scale_penultimate(sigxi_lambda_initial = params_initial_grid[j,], c = c, x = x, u = u, n_iter = n_iter_grid,
                                                      n_burn = n_burn_grid, sd_xi = sd_xi, sd_sigma = sd_sigma, min_lamda = min_lamda, max_lambda = max_lambda, prior_choice = prior)
      effectiveSize(t(params_samples))

    })
    mean(ess_start)
  })
  optim_sig_xi_sig <- grid[which.max(ess_grid),]

  # Run mcmc chains with the optimal variances obtained for different starting
  # points and asses convergence
  posterior_samples <- matrix(NA, nrow = 3, ncol = (n_iter - n_burn + 1) * n_starting_points) #preallocate space
  posterior_sigma <- list()
  posterior_xi <- list()
  posterior_lambda <- list()

  for(j in 1:n_starting_points) {
    posterior_samples[, ((n_iter - n_burn + 1) * (j - 1) + 1): ((n_iter - n_burn + 1) * j)] <- bayesian_estimation_gpd_scale_penultimate(sigxi_lambda_initial = params_initial_grid[j,], c = c, x = x, u = u, n_iter = n_iter,
                                                                                                                             n_burn = n_burn, sd_xi = optim_sig_xi_sig[[1]], sd_sigma = optim_sig_xi_sig[[2]], min_lamda = min_lamda, max_lambda = max_lambda, prior_choice = prior)
    posterior_sigma[[j]] <- mcmc(posterior_samples[1, ((n_iter - n_burn + 1) * (j - 1) + 1): ((n_iter - n_burn + 1) * j)])
    posterior_xi[[j]] <- mcmc(posterior_samples[2, ((n_iter - n_burn + 1) * (j - 1) + 1): ((n_iter - n_burn + 1) * j)])
    posterior_lambda[[j]] <- mcmc(posterior_samples[3, ((n_iter - n_burn + 1) * (j - 1) + 1): ((n_iter - n_burn + 1) * j)])
  }

  # Convergence diagnostics
  upper_ci_sigma <- gelman.diag(mcmc.list(posterior_sigma))[1]$psrf[2]
  #stopifnot('Chains with samples of sigma_u have not converge to the same distribution.' = upper_ci_sigma_u < 1.2)
  upper_ci_xi <- gelman.diag(mcmc.list(posterior_xi))[1]$psrf[2]
  #stopifnot('Chains with samples of xi have not converge to the same distribution.' = upper_ci_xi < 1.2)
  upper_ci_lambda <- gelman.diag(mcmc.list(posterior_lambda))[1]$psrf[2]

  if((upper_ci_sigma > 1.25) || (upper_ci_xi > 1.25) || (upper_ci_lambda > 1.25)){posterior_samples<- c(0,0)}

  return(posterior_samples)
}


## UNDER ASSUMPTION OF VALIDITY OF GPD FOR DIFFERENT LAMBDAS

# Import libraries
library(coda)
library(tidyr)

simulation_mcmc_scale <- function(xi, sigma,
                                  xi_init_min = -0.5, sigma_u_init_min = 0,
                                  xi_init_max = 1.5, sigma_u_init_max = 4,
                                  n_starting_points = 5, u,
                                  n_data = 500, n_iter = 1e4, n_burn = 1e3,
                                  prior = c('mdi', 'flat', 'jeffreys'),
                                  sd_xi = c(0.2, 0.15, 0.1, 0.05), sd_sigma = c(0.2, 0.15, 0.1, 0.05),
                                  min_lambda = -0.1, max_lambda = 3,
                                  scale_transformation = TRUE){

  prior <- match.arg(prior)

  # Generate data
  set.seed(123)
  x <- rgpd(n = n_data, scale = sigma, shape = xi, shift = u)
  if (scale_transformation == TRUE){
    # Data under scale transformation
    x <- log(x)
  }

  #Perform grid search to find optimal variances for the random walk
  n_iter_grid <- 1e3
  n_burn_grid <- 1e2
  grid <- crossing(sd_xi, sd_sigma)
  grid <- as.data.frame(grid)
  n_combinations <- dim(grid)[1]
  set.seed(111)
  params_initial_grid <- matrix(c(runif(n_starting_points, min = sigma_u_init_min, max = sigma_u_init_max),
                                  runif(n_starting_points, min = xi_init_min, max = xi_init_max),
                                  runif(n_starting_points, min = min_lambda, max = max_lambda)),
                                ncol = 3)

  ess_grid <- sapply(1:n_combinations, function(i){
    sd_xi <- grid[i, 1]
    sd_sigma <- grid[i, 2]
    ess_start <- sapply(1:n_starting_points, function(j){

      params_samples <- bayesian_estimation_gpd_scale(sigxi_lambda_initial = params_initial_grid[j,], x = x, u = u, n_iter = n_iter_grid,
                                                n_burn = n_burn_grid, sd_xi = sd_xi, sd_sigma = sd_sigma, min_lamda = min_lamda, max_lambda = max_lambda, prior_choice = prior)
      effectiveSize(t(params_samples))

    })
    mean(ess_start)
  })
  optim_sig_xi_sig <- grid[which.max(ess_grid),]

  # Run mcmc chains with the optimal variances obtained for different starting
  # points and asses convergence
  posterior_samples <- matrix(NA, nrow = 3, ncol = (n_iter - n_burn + 1) * n_starting_points) #preallocate space
  posterior_sigma <- list()
  posterior_xi <- list()
  posterior_lambda <- list()

  for(j in 1:n_starting_points) {
    posterior_samples[, ((n_iter - n_burn + 1) * (j - 1) + 1): ((n_iter - n_burn + 1) * j)] <- bayesian_estimation_gpd_scale(sigxi_lambda_initial = params_initial_grid[j,], x = x, u = u, n_iter = n_iter,
                                                                                                                       n_burn = n_burn, sd_xi = optim_sig_xi_sig[[1]], sd_sigma = optim_sig_xi_sig[[2]], min_lamda = min_lamda, max_lambda = max_lambda, prior_choice = prior)
    posterior_sigma[[j]] <- mcmc(posterior_samples[1, ((n_iter - n_burn + 1) * (j - 1) + 1): ((n_iter - n_burn + 1) * j)])
    posterior_xi[[j]] <- mcmc(posterior_samples[2, ((n_iter - n_burn + 1) * (j - 1) + 1): ((n_iter - n_burn + 1) * j)])
    posterior_lambda[[j]] <- mcmc(posterior_samples[3, ((n_iter - n_burn + 1) * (j - 1) + 1): ((n_iter - n_burn + 1) * j)])
  }

  # Convergence diagnostics
  upper_ci_sigma <- gelman.diag(mcmc.list(posterior_sigma))[1]$psrf[2]
  #stopifnot('Chains with samples of sigma_u have not converge to the same distribution.' = upper_ci_sigma_u < 1.2)
  upper_ci_xi <- gelman.diag(mcmc.list(posterior_xi))[1]$psrf[2]
  #stopifnot('Chains with samples of xi have not converge to the same distribution.' = upper_ci_xi < 1.2)
  upper_ci_lambda <- gelman.diag(mcmc.list(posterior_lambda))[1]$psrf[2]

  if((upper_ci_sigma > 1.25) || (upper_ci_xi > 1.25) || (upper_ci_lambda > 1.25)){posterior_samples<- c(0,0)}

  return(posterior_samples)
}
