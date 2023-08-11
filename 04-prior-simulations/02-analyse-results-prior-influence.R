############################
#SIMULATION 1
############################
library(ggplot2)

xi_sim <- -0.084
sigma_sim <- 0.48
u_sim <- 1.45
n_data_sim_vector <- c(45, 60, 80, 100, 125, 250, 500)
prior_sim_vector <- c('flat', 'mdi', 'jeffreys')

for (i in 1:length(n_data_sim_vector)){
  n_data_sim <- n_data_sim_vector[i]
  #Join posterior samples for different priors
  simulation_data_df <- data.frame()
  for (j in 1:length(prior_sim_vector)){
    prior_sim <- prior_sim_vector[j]

    file_name_sim <- paste('prior_influence',
                           prior_sim, 'prior',
                           n_data_sim, 'n_data', sep = '_')

    simulation_data <- read.csv(here::here('04-prior-simulations','outputs',file_name_sim),
                                colClasses=c("NULL", NA, NA), col.names = c('','sigma', 'xi'))

    simulation_data$prior <- rep(prior_sim, length(simulation_data$sigma))
    simulation_data_df <- rbind(simulation_data_df, simulation_data)
  }
  print(unique(simulation_data_df$prior))
  #Plot posterior xi-sigma
  plot_post <- ggplot(simulation_data_df, aes(x = xi, y = sigma, colour = prior)) +
    geom_density_2d()+
    xlab(expression(xi))+
    ylab(expression(sigma))+
    theme(#axis.title.x = element_blank(),
      #axis.title.y = element_blank(),
      axis.text.x = element_text(family = 'sans', size = 8),
      axis.text.y = element_text(family = 'sans', size = 8),
      legend.title = element_blank(),
      legend.text = element_text(family = 'sans', size = 8),
      legend.position = 'top',
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line.y = element_line(),
      axis.ticks.y = element_line())

  #save plot
  file_name <- paste('posterior-xi-sigma',
                     n_data_sim, 'n_data.png', sep = '-')

  ggsave(here::here('04-prior-simulations', 'figures', 'prior_influence', file_name),
         plot = plot_post,
         units = "px",
         width = 1197,
         height = 683,
         dpi = 300)
}

