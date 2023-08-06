###############
#SIMULATIONS
###############

# Import libraries
library(coda)
library(tidyr)

## PENULTIMATE APPROXIMATION
simulation_mcmc_scale_penultimate <- function(xi, sigma,
                                              xi_init_min = -0.5, sigma_init_min = 0,
                                              xi_init_max = 1, sigma_init_max = 4,
                                              n_starting_points = 5, u,
                                              n_data, n_iter = 1e4, n_burn = 1e3,
                                              prior = c('mdi', 'flat', 'jeffreys'),
                                              sd_xi = c(0.2, 0.15, 0.1, 0.05), sd_sigma = c(0.2, 0.15, 0.1, 0.05),
                                              min_lambda = 0, max_lambda = 3, step_lambda = 0.1,
                                              min_xi = -0.35, max_xi = 0.35, step_xi = 0.05,
                                              scale_transformation = TRUE){

  prior <- match.arg(prior)

  # Generate data
  set.seed(123)
  x <- rgpd(n = n_data, scale = sigma, shape = xi, shift = u)
  if (scale_transformation == TRUE){
    # Data under scale transformation
    x <- log(x)
    u <- log(u)
  }

  # Obtain estimate of c using weighted least squares regression
  # Initial values for MLE
  sigma_init_mle <- runif(1, sigma_init_min, sigma_init_max)
  xi_init_mle <- runif(1, min_xi, max_xi)
  lambda_init_mle <- runif(1, min_lambda, max_lambda)
  wlm_c <- estimate_c_penultimate(min_lambda = min_lambda, max_lambda = max_lambda, step_lambda = step_lambda,
                                  min_xi = min_xi, max_xi = max_xi, step_xi = step_xi, u = u, x = x, sigxi_lambda_init = c(sigma_init_mle, xi_init_mle, lambda_init_mle))
  c <- wlm_c[1]

  #Perform grid search to find optimal variances for the random walk
  n_iter_grid <- 1e3
  n_burn_grid <- 1e2
  grid <- crossing(sd_xi, sd_sigma)
  grid <- as.data.frame(grid)
  n_combinations <- dim(grid)[1]
  set.seed(111)
  params_initial_grid <- matrix(c(runif(n_starting_points, min = sigma_init_min, max = sigma_init_max),
                                  runif(n_starting_points, min = xi_init_min, max = xi_init_max),
                                  runif(n_starting_points, min = min_lambda, max = max_lambda)),
                                ncol = 3)

  ess_grid <- sapply(1:n_combinations, function(i){
    sd_xi <- grid[i, 1]
    sd_sigma <- grid[i, 2]
    ess_start <- sapply(1:n_starting_points, function(j){

      params_samples <- bayesian_estimation_gpd_scale_penultimate(sigxi_lambda_initial = params_initial_grid[j,], c = c, x = x, u = u, n_iter = n_iter_grid,
                                                      n_burn = n_burn_grid, sd_xi = sd_xi, sd_sigma = sd_sigma, min_lambda = min_lambda, max_lambda = max_lambda, prior_choice = prior)
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
                                                                                                                             n_burn = n_burn, sd_xi = optim_sig_xi_sig[[1]], sd_sigma = optim_sig_xi_sig[[2]], min_lambda = min_lambda, max_lambda = max_lambda, prior_choice = prior)
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

  #if((upper_ci_sigma > 1.25) || (upper_ci_xi > 1.25) || (upper_ci_lambda > 1.25)){posterior_samples<- c(0,0)}
  print(c(upper_ci_sigma, upper_ci_xi, upper_ci_lambda))

  return(posterior_samples)
}


## PENULTIMATE APPROXIMATION QUADRATIC
simulation_mcmc_scale_penultimate_quadratic <- function(xi, sigma,
                                                        xi_init_min = -0.5, sigma_init_min = 0,
                                                        xi_init_max = 1, sigma_init_max = 4,
                                                        n_starting_points = 5, u,
                                                        n_data, n_iter = 1e4, n_burn = 1e3,
                                                        prior = c('mdi', 'flat', 'jeffreys'),
                                                        sd_xi = c(0.2, 0.15, 0.1, 0.05), sd_sigma = c(0.2, 0.15, 0.1, 0.05),
                                                        min_lambda = 0, max_lambda = 3, step_lambda = 0.1,
                                                        min_xi = -0.35, max_xi = 0.35, step_xi = 0.05,
                                                        scale_transformation = TRUE){

  prior <- match.arg(prior)

  # Generate data
  set.seed(123)
  x <- rgpd(n = n_data, scale = sigma, shape = xi, shift = u)
  if (scale_transformation == TRUE){
    # Data under scale transformation
    x <- log(x)
    u <- log(u)
  }

  # Obtain estimate of c using weighted least squares regression
  # Initial values for MLE
  sigma_init_mle <- runif(1, sigma_init_min, sigma_init_max)
  xi_init_mle <- runif(1, min_xi, max_xi)
  lambda_init_mle <- runif(1, min_lambda, max_lambda)
  wlm_c <- estimate_c_penultimate_quadratic(min_lambda = min_lambda, max_lambda = max_lambda, step_lambda = step_lambda,
                                            min_xi = min_xi, max_xi = max_xi, step_xi = step_xi, u = u, x = x, sigxi_lambda_init = c(sigma_init_mle, xi_init_mle, lambda_init_mle))
  c1 <- wlm_c[1]
  c2 <- wlm_c[3]

  #Perform grid search to find optimal variances for the random walk
  n_iter_grid <- 1e3
  n_burn_grid <- 1e2
  grid <- crossing(sd_xi, sd_sigma)
  grid <- as.data.frame(grid)
  n_combinations <- dim(grid)[1]
  set.seed(111)
  params_initial_grid <- matrix(c(runif(n_starting_points, min = sigma_init_min, max = sigma_init_max),
                                  runif(n_starting_points, min = xi_init_min, max = xi_init_max),
                                  runif(n_starting_points, min = min_lambda, max = max_lambda)),
                                ncol = 3)

  ess_grid <- sapply(1:n_combinations, function(i){
    sd_xi <- grid[i, 1]
    sd_sigma <- grid[i, 2]
    ess_start <- sapply(1:n_starting_points, function(j){

      params_samples <- bayesian_estimation_gpd_scale_penultimate_quadratic(sigxi_lambda_initial = params_initial_grid[j,], c1 = c1, c2 = c2, x = x, u = u, n_iter = n_iter_grid,
                                                                            n_burn = n_burn_grid, sd_xi = sd_xi, sd_sigma = sd_sigma, min_lambda = min_lambda, max_lambda = max_lambda, prior_choice = prior)
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
    posterior_samples[, ((n_iter - n_burn + 1) * (j - 1) + 1): ((n_iter - n_burn + 1) * j)] <- bayesian_estimation_gpd_scale_penultimate_quadratic(sigxi_lambda_initial = params_initial_grid[j,], c1 = c1, c2 = c2, x = x, u = u, n_iter = n_iter,
                                                                                                                                         n_burn = n_burn, sd_xi = optim_sig_xi_sig[[1]], sd_sigma = optim_sig_xi_sig[[2]], min_lambda = min_lambda, max_lambda = max_lambda, prior_choice = prior)
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

  #if((upper_ci_sigma > 1.25) || (upper_ci_xi > 1.25) || (upper_ci_lambda > 1.25)){posterior_samples<- c(0,0)}
  print(c(upper_ci_sigma, upper_ci_xi, upper_ci_lambda))

  return(posterior_samples)
}


## UNDER ASSUMPTION OF VALIDITY OF GPD FOR DIFFERENT LAMBDAS

simulation_mcmc_scale <- function(xi, sigma,
                                  xi_init_min = -0.5, sigma_init_min = 0,
                                  xi_init_max = 1, sigma_init_max = 4,
                                  n_starting_points = 5, u,
                                  n_data, n_iter = 1e4, n_burn = 1e3,
                                  prior = c('mdi', 'flat', 'jeffreys'),
                                  sd_xi = c(0.2, 0.15, 0.1, 0.05), sd_sigma = c(0.2, 0.15, 0.1, 0.05),
                                  min_lambda = 0, max_lambda = 3,
                                  scale_transformation = TRUE){

  prior <- match.arg(prior)

  # Generate data
  set.seed(123)
  x <- rgpd(n = n_data, scale = sigma, shape = xi, shift = u)
  if (scale_transformation == TRUE){
    # Data under scale transformation
    x <- log(x)
    u <- log(u)
  }

  #Perform grid search to find optimal variances for the random walk
  n_iter_grid <- 1e3
  n_burn_grid <- 1e2
  grid <- crossing(sd_xi, sd_sigma)
  grid <- as.data.frame(grid)
  n_combinations <- dim(grid)[1]
  set.seed(111)
  params_initial_grid <- matrix(c(runif(n_starting_points, min = sigma_init_min, max = sigma_init_max),
                                  runif(n_starting_points, min = xi_init_min, max = xi_init_max),
                                  runif(n_starting_points, min = min_lambda, max = max_lambda)),
                                ncol = 3)

  ess_grid <- sapply(1:n_combinations, function(i){
    sd_xi <- grid[i, 1]
    sd_sigma <- grid[i, 2]
    ess_start <- sapply(1:n_starting_points, function(j){

      params_samples <- bayesian_estimation_gpd_scale(sigxi_lambda_initial = params_initial_grid[j,], x = x, u = u, n_iter = n_iter_grid,
                                                n_burn = n_burn_grid, sd_xi = sd_xi, sd_sigma = sd_sigma, min_lambda = min_lambda, max_lambda = max_lambda, prior_choice = prior)
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
                                                                                                                       n_burn = n_burn, sd_xi = optim_sig_xi_sig[[1]], sd_sigma = optim_sig_xi_sig[[2]], min_lambda = min_lambda, max_lambda = max_lambda, prior_choice = prior)
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

  #if((upper_ci_sigma > 1.25) || (upper_ci_xi > 1.25) || (upper_ci_lambda > 1.25)){posterior_samples<- c(0,0)}
  print(c(upper_ci_sigma, upper_ci_xi, upper_ci_lambda))

  return(posterior_samples)
}


# NOT TAKING INTO ACCOUNT SCALE UNCERTAINTY
simulation_mcmc <- function(xi, sigma,
                            xi_init_min = -0.5, sigma_init_min = 0,
                            xi_init_max = 1.5, sigma_init_max = 4,
                            n_starting_points = 5,
                            u,
                            n_data = 500, n_iter = 1e4, n_burn = 1e3,
                            prior = c('mdi', 'flat', 'jeffreys'),
                            sd_xi = c(0.2, 0.15, 0.1, 0.05), sd_sig = c(0.2, 0.15, 0.1, 0.05),
                            scale_transformation = TRUE){

  prior <- match.arg(prior)

  # Generate data
  set.seed(123)
  x <- rgpd(n = n_data, scale = sigma, shape = xi, shift = u)
  if (scale_transformation == TRUE){
    # Data under scale transformation
    x <- log(x)
    u <- log(u)
  }


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
  posterior_samples <- matrix(NA, nrow = 3, ncol = (n_iter - n_burn + 1) * n_starting_points) #preallocate space
  posterior_sigma <- list()
  posterior_xi <- list()



  for(j in 1:n_starting_points) {
    posterior_samples[1:2, ((n_iter - n_burn + 1) * (j - 1) + 1): ((n_iter - n_burn + 1) * j)] <- bayesian_estimation_gpd(sigxi_initial = params_initial_grid[j,], x = x, u = u, n_iter = n_iter,
                                                                                                                   n_burn = n_burn, sd_xi = optim_sig_xi_sig[[1]], sd_sig = optim_sig_xi_sig[[2]], prior_choice = prior)
    posterior_sigma[[j]] <- mcmc(posterior_samples[1, ((n_iter - n_burn + 1) * (j - 1) + 1): ((n_iter - n_burn + 1) * j)])
    posterior_xi[[j]] <- mcmc(posterior_samples[2, ((n_iter - n_burn + 1) * (j - 1) + 1): ((n_iter - n_burn + 1) * j)])
  }
  posterior_samples[3,] <- rep(1, (n_iter - n_burn + 1) * n_starting_points)

  # Convergence diagnostics
  upper_ci_sigma <- gelman.diag(mcmc.list(posterior_sigma))[1]$psrf[2]
  #stopifnot('Chains with samples of sigma have not converge to the same distribution.' = upper_ci_sigma < 1.2)
  upper_ci_xi <- gelman.diag(mcmc.list(posterior_xi))[1]$psrf[2]
  #stopifnot('Chains with samples of xi have not converge to the same distribution.' = upper_ci_xi < 1.2)

  #if((upper_ci_sigma_u > 1.25) || (upper_ci_xi > 1.25)){posterior_samples<- c(0,0)}
  print(c(upper_ci_sigma, upper_ci_xi))
  return(posterior_samples)

}


## PENULTIMATE APPROXIMATION UNCERTAINTY QUANTIFICATION
simulation_mcmc_scale_penultimate_uq <- function(xi, sigma,
                                              xi_init_min = -0.5, sigma_init_min = 0,
                                              xi_init_max = 1, sigma_init_max = 4,
                                              n_starting_points = 5, u,
                                              n_data, n_iter = 1e4, n_burn = 1e3,
                                              prior = c('mdi', 'flat', 'jeffreys'),
                                              sd_xi = c(0.2, 0.15, 0.1, 0.05), sd_sigma = c(0.2, 0.15, 0.1, 0.05),
                                              min_lambda = 0, max_lambda = 3, step_lambda = 0.1,
                                              min_xi = -0.35, max_xi = 0.35, step_xi = 0.05,
                                              scale_transformation = TRUE){

  prior <- match.arg(prior)

  # Generate data
  set.seed(123)
  x <- rgpd(n = n_data, scale = sigma, shape = xi, shift = u)
  if (scale_transformation == TRUE){
    # Data under scale transformation
    x <- log(x)
    u <- log(u)
  }

  # Obtain estimate of c using weighted least squares regression
  # Initial values for MLE
  sigma_init_mle <- runif(1, sigma_init_min, sigma_init_max)
  xi_init_mle <- runif(1, min_xi, max_xi)
  lambda_init_mle <- runif(1, min_lambda, max_lambda)
  wlm_c <- estimate_c_penultimate(min_lambda = min_lambda, max_lambda = max_lambda, step_lambda = step_lambda,
                                  min_xi = min_xi, max_xi = max_xi, step_xi = step_xi, u = u, x = x, sigxi_lambda_init = c(sigma_init_mle, xi_init_mle, lambda_init_mle))
  c <- wlm_c[1]
  sig_c <- wlm_c[2]
  print(c(c, sig_c))

  c_values <- seq(c-2*sig_c, c+2*sig_c, length.out = 7)
  #c_values <- c
  posterior_samples_join <- c()

  for(i in 1:length(c_values)){
    c <- c_values[i]
    #Perform grid search to find optimal variances for the random walk
    n_iter_grid <- 1e3
    n_burn_grid <- 1e2
    grid <- crossing(sd_xi, sd_sigma)
    grid <- as.data.frame(grid)
    n_combinations <- dim(grid)[1]
    set.seed(111)
    params_initial_grid <- matrix(c(runif(n_starting_points, min = sigma_init_min, max = sigma_init_max),
                                    runif(n_starting_points, min = xi_init_min, max = xi_init_max),
                                    runif(n_starting_points, min = min_lambda, max = max_lambda)),
                                  ncol = 3)

    ess_grid <- sapply(1:n_combinations, function(i){
      sd_xi <- grid[i, 1]
      sd_sigma <- grid[i, 2]
      ess_start <- sapply(1:n_starting_points, function(j){

        params_samples <- bayesian_estimation_gpd_scale_penultimate(sigxi_lambda_initial = params_initial_grid[j,], c = c, x = x, u = u, n_iter = n_iter_grid,
                                                                    n_burn = n_burn_grid, sd_xi = sd_xi, sd_sigma = sd_sigma, min_lambda = min_lambda, max_lambda = max_lambda, prior_choice = prior)
        effectiveSize(t(params_samples))

      })
      mean(ess_start)
    })
    optim_sig_xi_sig <- grid[which.max(ess_grid),]

    # Run mcmc chains with the optimal variances obtained for different starting
    # points and asses convergence
    posterior_samples <- matrix(NA, nrow = 4, ncol = (n_iter - n_burn + 1) * n_starting_points) #preallocate space
    posterior_sigma <- list()
    posterior_xi <- list()
    posterior_lambda <- list()

    for(j in 1:n_starting_points) {
      posterior_samples[1:3, ((n_iter - n_burn + 1) * (j - 1) + 1): ((n_iter - n_burn + 1) * j)] <- bayesian_estimation_gpd_scale_penultimate(sigxi_lambda_initial = params_initial_grid[j,], c = c, x = x, u = u, n_iter = n_iter,
                                                                                                                                           n_burn = n_burn, sd_xi = optim_sig_xi_sig[[1]], sd_sigma = optim_sig_xi_sig[[2]], min_lambda = min_lambda, max_lambda = max_lambda, prior_choice = prior)
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

    #if((upper_ci_sigma > 1.25) || (upper_ci_xi > 1.25) || (upper_ci_lambda > 1.25)){posterior_samples<- c(0,0)}
    print(c(upper_ci_sigma, upper_ci_xi, upper_ci_lambda))

    posterior_samples[4,] <- rep(paste(c, sep = ''), (n_iter - n_burn + 1) * n_starting_points)
    posterior_samples_join <- cbind(posterior_samples_join, posterior_samples)
  }
  return(posterior_samples_join)
}

## PENULTIMATE APPROXIMATION UNCERTAINTY QUANTIFICATION
simulation_mcmc_scale_penultimate_uq_c_estimates <- function(xi, sigma,
                                                 xi_init_min = -0.5, sigma_init_min = 0,
                                                 xi_init_max = 1, sigma_init_max = 4,
                                                 n_starting_points = 5, u,
                                                 n_data, n_iter = 1e4, n_burn = 1e3,
                                                 prior = c('mdi', 'flat', 'jeffreys'),
                                                 sd_xi = c(0.2, 0.15, 0.1, 0.05), sd_sigma = c(0.2, 0.15, 0.1, 0.05),
                                                 min_lambda = 0, max_lambda = 3, step_lambda = 0.1,
                                                 min_xi = -0.35, max_xi = 0.35, step_xi = 0.05,
                                                 scale_transformation = TRUE){

  prior <- match.arg(prior)

  # Generate data
  set.seed(123)
  x <- rgpd(n = n_data, scale = sigma, shape = xi, shift = u)
  if (scale_transformation == TRUE){
    # Data under scale transformation
    x <- log(x)
    u <- log(u)
  }

  # Obtain estimate of c using weighted least squares regression
  # Initial values for MLE
  sigma_init_mle <- runif(1, sigma_init_min, sigma_init_max)
  xi_init_mle <- runif(1, min_xi, max_xi)
  lambda_init_mle <- runif(1, min_lambda, max_lambda)
  wlm_c <- estimate_c_penultimate(min_lambda = min_lambda, max_lambda = max_lambda, step_lambda = step_lambda,
                                  min_xi = min_xi, max_xi = max_xi, step_xi = step_xi, u = u, x = x, sigxi_lambda_init = c(sigma_init_mle, xi_init_mle, lambda_init_mle))
  c <- wlm_c[1]
  sig_c <- wlm_c[2]
  return(c(c, sig_c))
}
