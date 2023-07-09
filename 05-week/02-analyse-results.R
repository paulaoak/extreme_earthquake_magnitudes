###################################
#READ DATA AND COMPUTE QUANTILES
###################################


library(ggplot2)
############################
#SIMULATION 1 FLAT
############################
xi_sim <- -0.084
sigma_u_sim <- 0.48
u_sim <- 1
n_data_sim <- 500
v_sim <- rep(c(1.1, 1.45), each = n_data_sim/2)

prior_sim_vector <- c('flat', 'jeffreys', 'mdi')

for (i in 1:length(prior_sim_vector)){
  prior_sim <- prior_sim_vector[i]
  file_name_sim <- paste('simulation', prior_sim, 'prior',
                         'xi', xi_sim,
                         'sigma_u', sigma_u_sim,
                         'u', u_sim,
                         'v_min', min(v_sim),
                         'v_max', max(v_sim),
                         n_data_sim, 'n_data', sep = '_')

  simulation_data_flat_1 <- read.csv(here::here('04-week','simulation_mcmc', 'outputs_variable_new',file_name_sim),
                                     colClasses=c("NULL", NA, NA, NA, NA), col.names = c('','sigma_all', 'xi_all', 'sigma_high_thres', 'xi_high_thres'))

  simulation_1_data <- data.frame(sigma = c(simulation_data_flat_1$sigma_all, simulation_data_flat_1$sigma_high_thres),
                                  xi = c(simulation_data_flat_1$xi_all, simulation_data_flat_1$xi_high_thres),
                                  threshold = rep(c('Variable', 'Constant'), each = length(simulation_data_flat_1$sigma_all)))

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
