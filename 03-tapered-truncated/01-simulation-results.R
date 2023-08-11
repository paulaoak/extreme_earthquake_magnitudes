############################
#SIMULATION 1
############################
xi_sim <- -0.084
sigma_sim <- 0.48
u_sim <- 1.45
#n_data_sim_vector <- c(60, 125, 250, 500)
#n_data_sim_vector <- c(60, 75, 125, 250, 500)
n_data_sim_vector <- c(75)

for(i in 1:length(n_data_sim_vector)){
  n_data_sim <- n_data_sim_vector[i]


  #Simulation scale under assumption of validity of GPD model for different
  simulation_results_2 <- simulation_mcmc_tapered(xi = xi_sim, sigma = sigma_sim,
                                                beta_init_min = 1, theta_init_min = 1,
                                                beta_init_max = 5, theta_init_max = 3,
                                                n_starting_points = 5, u = u_sim,
                                                n_data = n_data_sim, n_iter = 1e4, n_burn = 1e3,
                                                sd_theta = c(0.1), sd_beta = c(0.1))

  file_name_sim_2 <- paste('tapered_simulation',
                           n_data_sim, 'n_data',
                           sep = '_')


  write.csv(t(simulation_results_2),
            here::here('03-tapered-truncated', 'outputs', file_name_sim_2))

  #Simulation scale penultimate
  simulation_results <- simulation_mcmc_truncated_gr(xi = xi_sim, sigma = sigma_sim,
                                                     n_starting_points = 5, u = u_sim,
                                                     n_data = n_data_sim, n_iter = 1e4, n_burn = 1e3)

  file_name_sim <- paste('truncated_gr_simulation',
                         n_data_sim, 'n_data',
                         sep = '_')


  write.csv(t(simulation_results),
            here::here('03-tapered-truncated', 'outputs',file_name_sim))

  #simulation without taking into account scale uncertainty
  simulation_results_3 <- simulation_mcmc_gpd(xi = xi_sim, sigma = sigma_sim,
                                          xi_init_min = -0.5, sigma_init_min = 0,
                                          xi_init_max = 1, sigma_init_max = 4,
                                          n_starting_points = 5,
                                          u = u_sim,
                                          n_data = 500, n_iter = 1e4, n_burn = 1e3,
                                          sd_xi = c(0.2, 0.15, 0.1, 0.05), sd_sig = c(0.2, 0.15, 0.1, 0.05))

  file_name_sim_3 <- paste('gpd_simulation', n_data_sim, 'n_data', sep = '_')


  write.csv(t(simulation_results_3),
            here::here('03-tapered-truncated', 'outputs', file_name_sim_3))

}
