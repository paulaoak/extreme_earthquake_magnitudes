###############
#SIMULATIONS
###############

# Import libraries
library(coda)
library(tidyr)

simulation_mcmc_1 <- function(mmax, mean,
                            n_starting_points = 5,
                            u, n_data = 500, n_iter = 1e4, n_burn = 1e3,
                            prior = c('flat-flat', 'flat-gamma'),
                            b_value, epsilon, upper_mmax = NULL, alpha = NULL, beta = NULL,
                            sd_mmax = c(0.1, 0.05), sd_mean = c(0.9)){

  prior <- match.arg(prior)

  # Generate data
  set.seed(123)
  x <- rgpd_mmax_mean(n = n_data, mmax = mmax, mean = mean, shift = u)

  #Perform grid search to find optimal variances for the random walk
  n_iter_grid <- 1e3
  n_burn_grid <- 1e2
  grid <- crossing(sd_mmax, sd_mean) #grid of possible variances
  grid <- as.data.frame(grid)
  n_combinations <- dim(grid)[1]

  #Obtain grid for starting parameters
  if(prior == 'flat-flat'){
    set.seed(111)
    params_initial_grid <- matrix(c(runif(n_starting_points, min = b_value + epsilon, max = upper_mmax),
                                    runif(n_starting_points, min = u, max = b_value + epsilon)),
                                  ncol = 2)
  }
  else{
    set.seed(111)
    params_initial_grid <- matrix(c(runif(n_starting_points, min = b_value + epsilon, max = 4* b_value + epsilon),
                                    runif(n_starting_points, min = u, max = b_value + epsilon)),
                                  ncol = 2)
  }

  ess_grid <- sapply(1:n_combinations, function(i){
    sd_mmax <- grid[i, 1]
    sd_mean <- grid[i, 2]
    ess_start <- sapply(1:n_starting_points, function(j){
      params_samples <- bayesian_estimation_gpd_mmax_mean(mmax_mean_initial = params_initial_grid[j,], x = x, u = u, n_iter = n_iter_grid,
                                                n_burn = n_burn_grid, sd_mmax = sd_mmax, sd_mean = sd_mean, prior_choice = prior,
                                                b_value = b_value, epsilon = epsilon, alpha = alpha, beta = beta, upper_mmax = upper_mmax)
      effectiveSize(t(params_samples))


    })
    mean(ess_start)
  })
  optim_sig_mmax_sig_mean <- grid[which.max(ess_grid),]

  # Run mcmc chains with the optimal variances obtained for different starting
  # points and asses convergence
  posterior_samples <- matrix(NA, nrow = 2, ncol = (n_iter - n_burn + 1) * n_starting_points) #preallocate space
  posterior_mmax <- list()
  posterior_mean <- list()


  for(j in 1:n_starting_points) {
    posterior_samples[, ((n_iter - n_burn + 1) * (j - 1) + 1): ((n_iter - n_burn + 1) * j)] <- bayesian_estimation_gpd_mmax_mean(mmax_mean_initial = params_initial_grid[j,], x = x, u = u, n_iter = n_iter,
                                                                                                                       n_burn = n_burn, sd_mmax = optim_sig_mmax_sig_mean[[1]], sd_mean = optim_sig_mmax_sig_mean[[2]], prior_choice = prior,
                                                                                                                       b_value = b_value, epsilon = epsilon, alpha = alpha, beta = beta, upper_mmax = upper_mmax)
    posterior_mmax[[j]] <- mcmc(posterior_samples[1, ((n_iter - n_burn + 1) * (j - 1) + 1): ((n_iter - n_burn + 1) * j)])
    posterior_mean[[j]] <- mcmc(posterior_samples[2, ((n_iter - n_burn + 1) * (j - 1) + 1): ((n_iter - n_burn + 1) * j)])
  }

  traceplot(posterior_mmax)

  # Convergence diagnostics
  upper_ci_mmax <- gelman.diag(mcmc.list(posterior_mmax))[1]$psrf[2]
  #stopifnot('Chains with samples of sigma_u have not converge to the same distribution.' = upper_ci_sigma_u < 1.2)
  upper_ci_mean <- gelman.diag(mcmc.list(posterior_mean))[1]$psrf[2]
  #stopifnot('Chains with samples of xi have not converge to the same distribution.' = upper_ci_xi < 1.2)

  if((upper_ci_mmax > 1.25) || (upper_ci_mean > 1.25)){posterior_samples<- c(0,0)}
  print(c(upper_ci_mmax, upper_ci_mean))
  return(posterior_samples)


}


############################
#SIMULATION 1
############################
mmax_sim <- 7.6
mean_sim <- 2.1
u_sim <- 1.45
b_value_sim <- 1.8
epsilon_sim <- 1.5
upper_mmax_sim <- 8.5
alpha_sim <- 1
beta_sim <- 0.5
prior_sim_vector <- c('flat-gamma', 'flat-flat')
n_data_sim_vector <- c(60, 125, 250, 500)


for(i in 1:length(n_data_sim_vector)){
  n_data_sim <- n_data_sim_vector[i]
  for(j in 1:length(prior_sim_vector)){
    prior_sim <- prior_sim_vector[j]

    simulation_results <- simulation_mcmc_1(mmax = mmax_sim, mean = mean_sim,
                                            u = u_sim, n_data = 500, n_iter = 1e4, n_burn = 1e3,
                                            prior = prior_sim,
                                            b_value = b_value_sim, epsilon = epsilon_sim, upper_mmax = upper_mmax_sim, alpha = alpha_sim, beta = beta_sim)


    file_name_sim <- paste('mmax', mmax_sim,
                           'mean', mean_sim,
                           'u', u_sim,
                           prior_sim, 'prior',
                           n_data_sim, 'n_data', sep = '_')

    write.csv(t(simulation_results),
              here::here('05-week','outputs',file_name_sim))
  }
}
