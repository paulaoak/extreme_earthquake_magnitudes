############################
#SIMULATION 1
############################
xi_sim <- -0.084
sigma_sim <- 0.48
u_sim <- 1.45
n_data_sim_vector <- c(45, 60, 80, 100, 125, 250, 500)
prior_sim_vector <- c('flat', 'mdi', 'jeffreys')

for (i in 1:length(n_data_sim_vector)){
  n_data_sim <- n_data_sim_vector[i]
  for(j in 1:length(prior_sim_vector)){
    prior_sim <- prior_sim_vector[j]

    simulation_results <- simulation_mcmc(xi = xi_sim, sigma = sigma_sim,
                                          xi_init_min = -0.5, sigma_init_min = 0,
                                          xi_init_max = 1.5, sigma_init_max = 4,
                                          n_starting_points = 5,
                                          u = u_sim,
                                          n_data = n_data_sim, n_iter = 1e4, n_burn = 1e3,
                                          prior = prior_sim)

    file_name_sim <- paste('prior_influence',
                             prior_sim, 'prior',
                             n_data_sim, 'n_data', sep = '_')

    write.csv(t(simulation_results),
              here::here('04-prior-simulations', 'outputs',file_name_sim))

  }
}
