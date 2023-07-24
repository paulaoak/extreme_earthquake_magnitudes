############################
#SIMULATION 1
############################
xi_sim <- -0.084
sigma_sim <- 0.48
u_sim <- 1.45
#n_data_sim_vector <- c(60, 125, 250, 500)
n_data_sim_vector <- c(125)
prior_sim_vector <- c('flat', 'mdi', 'jeffreys')
scale_transformation_sim_vector <- c(TRUE, FALSE)
#scale_transformation_sim_vector <- c(FALSE)

for(i in 1:length(n_data_sim_vector)){
  n_data_sim <- n_data_sim_vector[i]
  for(j in 1:length(prior_sim_vector)){
    prior_sim <- prior_sim_vector[j]
    for(k in 1:length(scale_transformation_sim_vector)){
      scale_transformation_sim <- scale_transformation_sim_vector[k]
      #Simulation scale penultimate
      simulation_results <- simulation_mcmc_scale_penultimate(xi = xi_sim, sigma = sigma_sim,
                                                              xi_init_min = -0.5, sigma_init_min = 0,
                                                              xi_init_max = 1, sigma_init_max = 4,
                                                              n_starting_points = 5, u = u_sim,
                                                              n_data = n_data_sim, n_iter = 1e4, n_burn = 1e3,
                                                              prior = prior_sim,
                                                              sd_xi = c(0.2, 0.15, 0.1, 0.05), sd_sigma = c(0.2, 0.15, 0.1, 0.05),
                                                              min_lambda = -0.1, max_lambda = 2, step_lambda = 0.1,
                                                              min_xi = -0.35, max_xi = 0.35, step_xi = 0.05,
                                                              scale_transformation = scale_transformation_sim)

      file_name_sim <- paste('simulation_penultimate', 'xi', xi_sim,
                             'sigma', sigma_sim, 'u', u_sim,
                             prior_sim, 'prior', n_data_sim, 'n_data',
                             'transform_data',scale_transformation_sim, sep = '_')


      write.csv(t(simulation_results),
                here::here('06-measurement-scale-uncertainty', 'outputs',file_name_sim))

      #Simulation scale under assumption of validity of GPD model for different
      simulation_results_2 <- simulation_mcmc_scale(xi = xi_sim, sigma = sigma_sim,
                                                    xi_init_min = -0.5, sigma_init_min = 0,
                                                    xi_init_max = 1, sigma_init_max = 4,
                                                    n_starting_points = 5, u = u_sim,
                                                    n_data = n_data_sim, n_iter = 1e4, n_burn = 1e3,
                                                    prior = prior_sim,
                                                    sd_xi = c(0.2, 0.15, 0.1, 0.05), sd_sigma = c(0.2, 0.15, 0.1, 0.05),
                                                    min_lambda = -0.1, max_lambda = 3,
                                                    scale_transformation = scale_transformation_sim)

      file_name_sim_2 <- paste('simulation', 'xi', xi_sim,
                             'sigma', sigma_sim, 'u', u_sim,
                             prior_sim, 'prior', n_data_sim, 'n_data',
                             'transform_data',scale_transformation_sim, sep = '_')


      write.csv(t(simulation_results_2),
                here::here('06-measurement-scale-uncertainty', 'outputs', file_name_sim_2))
      }
  }
}

