###################################
#READ DATA AND COMPUTE QUANTILES
###################################


library(ggplot2)

############################
#SIMULATION
############################
mmax_sim <- 4.6
mean_sim <- 2.1
u_sim <- 1.45
b_value_sim <- 1.8
epsilon_sim <- 0.8
upper_mmax_sim <- 6
prior_sim_vector <- c('flat-flat', 'flat-gamma')
n_data_sim_vector <- c(75, 125, 250, 500)


for(i in 1:length(n_data_sim_vector)){
  n_data_sim <- n_data_sim_vector[i]
  for(j in 1:length(prior_sim_vector)){
    prior_sim <- prior_sim_vector[j]

    simulation_results <- simulation_mcmc_1(mmax = mmax_sim, mean = mean_sim,
                                            u = u_sim, n_data = 500, n_iter = 1e4, n_burn = 1e3,
                                            prior = prior_sim,
                                            b_value = b_value_sim, epsilon = epsilon_sim, upper_mmax = upper_mmax_sim)


    file_name_sim <- paste('mmax', mmax_sim,
                           'mean', mean_sim,
                           'u', u_sim,
                           prior_sim, 'prior',
                           n_data_sim, 'n_data', sep = '_')

    simulation_data <- read.csv(here::here('05-week','outputs',file_name_sim),
                                colClasses=c("NULL", NA, NA), col.names = c('','mmax', 'mean'))

    simulation_data_df <- data.frame(mmax = c(simulation_data$mmax, simulation_data_flat_1$sigma_high_thres),
                                    xi = c(simulation_data_flat_1$xi_all, simulation_data_flat_1$xi_high_thres),
                                    threshold = rep(c('Variable', 'Constant'), each = length(simulation_data_flat_1$sigma_all)))

  }
}

  if(prior_sim == 'flat'){prior_sim = 'Flat'}
  if(prior_sim == 'jeffreys'){prior_sim = 'Jeffreys'}
  if(prior_sim == 'mdi'){prior_sim = 'MDI'}
  plot_post <- ggplot(simulation_1_data, aes(x = xi, y = sigma, colour = threshold)) +
    geom_density_2d()+
    xlab(expression(xi))+
    ylab(expression(sigma))+
    ggtitle(paste(prior_sim, 'prior', sep = ' '))+
    theme(#axis.title.x = element_blank(),
      #axis.title.y = element_blank(),
      axis.text.x = element_text(family = 'sans'),
      axis.text.y = element_text(family = 'sans'),
      legend.title = element_blank(),
      legend.text = element_text(family = 'sans', size = 7),
      legend.position = 'top',
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line.y = element_line(),
      axis.ticks.y = element_line())

  file_name <- paste('posterior-xi-sigma-title',
                     'sigma', sigma_u_sim,
                     'xi', xi_sim,
                     prior_sim, 'prior.png', sep = '-')

  ggsave(here::here('04-week', 'simulation_mcmc', 'figures','posterior_xi_sigma_threshold_poster', file_name),
         plot = plot_post,
         units = "px",
         width = 1197,
         height = 683,
         dpi = 300)
}
