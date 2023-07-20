###################################
#READ DATA AND COMPUTE QUANTILES
###################################


library(ggplot2)

############################
#SIMULATION
############################
mmax_sim <- 7.6
mean_sim <- 2.1
u_sim <- 1.45
b_value_sim <- 1.8
epsilon_sim <- 1.5
upper_mmax_sim <- 8.5
alpha_sim <- 1
beta_sim <- 0.5
prior_sim <- 'flat-flat'
n_data_sim_vector <- c(60, 125, 250, 500)

#obtain data frame with results from different simulations
simulation_data_df <- data.frame()
for (i in 1:length(n_data_sim_vector)){
  n_data_sim <- n_data_sim_vector[i]
  file_name_sim <- paste('mmax', mmax_sim,
                         'mean', mean_sim,
                         'u', u_sim,
                         prior_sim, 'prior',
                         n_data_sim, 'n_data', sep = '_')

  simulation_data <- read.csv(here::here('05-week','outputs_2',file_name_sim),
                              colClasses=c("NULL", NA, NA), col.names = c('','mmax', 'mean'))

  simulation_data$type <- rep(paste(n_data_sim, 'obs', sep = ' '), length(simulation_data$mean))

  simulation_data_df <- rbind(simulation_data_df, simulation_data)

}

simulation_data_df$type <- factor(simulation_data_df$type, levels = c('500 obs', '250 obs', '125 obs','60 obs'))
#if(prior_sim == 'flat'){prior_sim = 'Flat'}
#if(prior_sim == 'jeffreys'){prior_sim = 'Jeffreys'}
#if(prior_sim == 'mdi'){prior_sim = 'MDI'}
plot_post <- ggplot(simulation_data_df, aes(x = mean, y = mmax, colour = type)) +
  geom_density_2d()+
  xlab('mean')+
  ylab(expression(M_max))+
  ggtitle(paste('Posterior density for', prior_sim, 'prior', sep = ' '))+
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

file_name <- paste('posterior-mmax-mean',
                   'mmax', mmax_sim,
                   'mean', mean_sim,
                   prior_sim, 'prior.png', sep = '-')

ggsave(here::here('05-week','figures','posterior_mmax_mean_different_obs', file_name),
       plot = plot_post,
       units = "px",
       width = 1197,
       height = 683,
       dpi = 300)

