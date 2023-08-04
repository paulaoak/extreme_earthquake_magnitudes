###################################
#READ DATA AND COMPUTE QUANTILES
###################################


library(ggplot2)
library(urbnthemes)

#function to compute quantile from mean and mmax
quantile_posterior_calculation <- function(mmax, mean, quantile, u){
  p <- 1 - quantile
  mmax = mmax - u
  mean = mean - u
  xi = mean / (mean - mmax)
  sigma = mean * mmax / (mmax - mean)

  quantiles_posterior <- u + sigma/xi * (p^(-xi) - 1)
}

############################
#Analysis and visualization 1
############################
xi_sim <- -0.084
sigma_sim <- 0.48
u_sim <- 1.45
n_data_sim_vector <- c(60, 75, 125, 250, 500)
prior_sim_vector <- c('flat')
scale_transformation_sim_vector <- c(TRUE, FALSE)

#obtain data frame with results from different simulations
simulation_data_df <- data.frame()
for (j in 1:length(prior_sim_vector)){
  prior_sim <- prior_sim_vector[j]
  for (i in 1:length(n_data_sim_vector)){
    n_data_sim <- n_data_sim_vector[i]

    file_name_sim_penultimate <- paste('simulation_penultimate', 'xi', xi_sim,
                                       'sigma', sigma_sim, 'u', u_sim,
                                       prior_sim, 'prior', n_data_sim, 'n_data',
                                       'transform_data',scale_transformation_sim, sep = '_')

    file_name_sim_2 <- paste('simulation', 'xi', xi_sim,
                             'sigma', sigma_sim, 'u', u_sim,
                             prior_sim, 'prior', n_data_sim, 'n_data',
                             'transform_data',scale_transformation_sim, sep = '_')

    simulation_data <- read.csv(here::here('06-measurement-scale-uncertainty','outputs', 'figures',file_name_sim),
                                colClasses=c("NULL", NA, NA, NA), col.names = c('','sigma', 'xi', 'lambda'))

    #add a column to indicate the number of observations used for the simulation
    simulation_data$type <- rep(paste(n_data_sim, 'obs', sep = ' '), length(simulation_data$mean))

    #add a column to compute the 0.5, 0.75, 0.9 and 0.95 quantiles
    simulation_data$quantile_0.5 <- quantile_posterior_calculation(mmax = simulation_data$mmax, mean = simulation_data$mean, u = u_sim, quantile = 0.5)
    simulation_data$quantile_0.75 <- quantile_posterior_calculation(mmax = simulation_data$mmax, mean = simulation_data$mean, u = u_sim, quantile = 0.75)
    simulation_data$quantile_0.9 <- quantile_posterior_calculation(mmax = simulation_data$mmax, mean = simulation_data$mean, u = u_sim, quantile = 0.9)
    simulation_data$quantile_0.95 <- quantile_posterior_calculation(mmax = simulation_data$mmax, mean = simulation_data$mean, u = u_sim, quantile = 0.95)

    simulation_data_df <- rbind(simulation_data_df, simulation_data)

  }

}


simulation_data_df$type <- factor(simulation_data_df$type, levels = c('500 obs', '250 obs', '125 obs','60 obs'))
#if(prior_sim == 'flat'){prior_sim = 'Flat'}
#if(prior_sim == 'jeffreys'){prior_sim = 'Jeffreys'}
#if(prior_sim == 'mdi'){prior_sim = 'MDI'}

#Plot posterior density
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

#Plot 0.5 quantile
quantile<- 0.5
true_quantile <- quantile_posterior_calculation(mmax = mmax_sim ,mean = mean_sim, quantile =quantile, u = u_sim)
plot_quantile <- ggplot(simulation_data_df, aes(x = quantile_0.5, colour = type, fill = type)) +
  geom_density(alpha = 0.4)+
  geom_segment(aes(x = true_quantile, y = 0, xend = true_quantile, yend = Inf, linetype = "True quantile"), color = 'black')+
  ggtitle(paste('Posterior distribution of quantile', quantile, sep = ' '))+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(family = 'sans'),
        axis.text.y = element_text(family = 'sans'),
        title = element_text(family = 'sans', size = 8),
        legend.title = element_blank(),
        legend.text = element_text(family = 'sans', size = 7),
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        axis.line.y = element_line(),
        axis.ticks.y = element_line(),
        axis.line.x = element_line(),
        axis.ticks.x = element_line())+
  #legend.direction='vertical')+
  #legend.position = 'right')+
  expand_limits(y = 0) +
  coord_cartesian(expand = FALSE, clip = "off")


#save plot
file_name <- paste('posterior_0.5_quantile',
                   prior_sim, 'prior',
                   mmax_sim, 'mmax',
                   mean_sim, 'mean.png', sep = '-')

ggsave(here::here('05-week', 'figures','posterior_mmax_mean_different_obs', file_name),
       plot = plot_quantile,
       units = "px",
       width = 1197,
       height = 683,
       dpi = 300)

#Plot 0.75 quantile
quantile<- 0.75
true_quantile <- quantile_posterior_calculation(mmax = mmax_sim ,mean = mean_sim, quantile =quantile, u = u_sim)
plot_quantile <- ggplot(simulation_data_df, aes(x = quantile_0.75, colour = type, fill = type)) +
  geom_density(alpha = 0.4)+
  geom_segment(aes(x = true_quantile, y = 0, xend = true_quantile, yend = Inf, linetype = "True quantile"), color = 'black')+
  ggtitle(paste('Posterior distribution of quantile', quantile, sep = ' '))+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(family = 'sans'),
        axis.text.y = element_text(family = 'sans'),
        title = element_text(family = 'sans', size = 8),
        legend.title = element_blank(),
        legend.text = element_text(family = 'sans', size = 7),
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        axis.line.y = element_line(),
        axis.ticks.y = element_line(),
        axis.line.x = element_line(),
        axis.ticks.x = element_line())+
  #legend.direction='vertical')+
  #legend.position = 'right')+
  expand_limits(y = 0) +
  coord_cartesian(expand = FALSE, clip = "off")


#save plot
file_name <- paste('posterior_0.75_quantile',
                   prior_sim, 'prior',
                   mmax_sim, 'mmax',
                   mean_sim, 'mean.png', sep = '-')

ggsave(here::here('05-week', 'figures','posterior_mmax_mean_different_obs', file_name),
       plot = plot_quantile,
       units = "px",
       width = 1197,
       height = 683,
       dpi = 300)


#Plot 0.9 quantile
quantile<- 0.9
true_quantile <- quantile_posterior_calculation(mmax = mmax_sim ,mean = mean_sim, quantile =quantile, u = u_sim)
plot_quantile <- ggplot(simulation_data_df, aes(x = quantile_0.9, colour = type, fill = type)) +
  geom_density(alpha = 0.4)+
  geom_segment(aes(x = true_quantile, y = 0, xend = true_quantile, yend = Inf, linetype = "True quantile"), color = 'black')+
  ggtitle(paste('Posterior distribution of quantile', quantile, sep = ' '))+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(family = 'sans'),
        axis.text.y = element_text(family = 'sans'),
        title = element_text(family = 'sans', size = 8),
        legend.title = element_blank(),
        legend.text = element_text(family = 'sans', size = 7),
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        axis.line.y = element_line(),
        axis.ticks.y = element_line(),
        axis.line.x = element_line(),
        axis.ticks.x = element_line())+
  #legend.direction='vertical')+
  #legend.position = 'right')+
  expand_limits(y = 0) +
  coord_cartesian(expand = FALSE, clip = "off")


#save plot
file_name <- paste('posterior_0.9_quantile',
                   prior_sim, 'prior',
                   mmax_sim, 'mmax',
                   mean_sim, 'mean.png', sep = '-')

ggsave(here::here('05-week', 'figures','posterior_mmax_mean_different_obs', file_name),
       plot = plot_quantile,
       units = "px",
       width = 1197,
       height = 683,
       dpi = 300)

#Plot 0.95 quantile
quantile<- 0.95
true_quantile <- quantile_posterior_calculation(mmax = mmax_sim ,mean = mean_sim, quantile =quantile, u = u_sim)
plot_quantile <- ggplot(simulation_data_df, aes(x = quantile_0.95, colour = type, fill = type)) +
  geom_density(alpha = 0.4)+
  geom_segment(aes(x = true_quantile, y = 0, xend = true_quantile, yend = Inf, linetype = "True quantile"), color = 'black')+
  ggtitle(paste('Posterior distribution of quantile', quantile, sep = ' '))+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(family = 'sans'),
        axis.text.y = element_text(family = 'sans'),
        title = element_text(family = 'sans', size = 8),
        legend.title = element_blank(),
        legend.text = element_text(family = 'sans', size = 7),
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        axis.line.y = element_line(),
        axis.ticks.y = element_line(),
        axis.line.x = element_line(),
        axis.ticks.x = element_line())+
  #legend.direction='vertical')+
  #legend.position = 'right')+
  expand_limits(y = 0) +
  coord_cartesian(expand = FALSE, clip = "off")


#save plot
file_name <- paste('posterior_0.95_quantile',
                   prior_sim, 'prior',
                   mmax_sim, 'mmax',
                   mean_sim, 'mean.png', sep = '-')

ggsave(here::here('05-week', 'figures','posterior_mmax_mean_different_obs', file_name),
       plot = plot_quantile,
       units = "px",
       width = 1197,
       height = 683,
       dpi = 300)


############################
#Analysis and visualization 2
############################
mmax_sim <- 7.6
mean_sim <- 2.1
u_sim <- 1.45
b_value_sim <- 1.8
epsilon_sim <- 1.5
upper_mmax_sim <- 8.5
alpha_sim_1 <- 1
beta_sim_1 <- 0.5
alpha_sim_2 <- 1
beta_sim_2 <- 0.5
prior_sim_vector <- c('unif-unif', 'unif-gamma', 'flat-gamma', 'gamma-gamma', 'flat-flat')
#n_data_sim_vector <- c(60, 125, 250, 500)
n_data_sim_vector <- c(25, 35, 45, 55)

for (j in 1:length(prior_sim_vector)){
  prior_sim <- prior_sim_vector[j]
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

    #add a column to indicate the number of observations used for the simulation
    simulation_data$type <- rep(paste(n_data_sim, 'obs', sep = ' '), length(simulation_data$mean))

    #add a column to compute the 0.5, 0.75, 0.9 and 0.95 quantiles
    simulation_data$quantile_0.5 <- quantile_posterior_calculation(mmax = simulation_data$mmax, mean = simulation_data$mean, u = u_sim, quantile = 0.5)
    simulation_data$quantile_0.75 <- quantile_posterior_calculation(mmax = simulation_data$mmax, mean = simulation_data$mean, u = u_sim, quantile = 0.75)
    simulation_data$quantile_0.9 <- quantile_posterior_calculation(mmax = simulation_data$mmax, mean = simulation_data$mean, u = u_sim, quantile = 0.9)
    simulation_data$quantile_0.95 <- quantile_posterior_calculation(mmax = simulation_data$mmax, mean = simulation_data$mean, u = u_sim, quantile = 0.95)

    simulation_data_df <- rbind(simulation_data_df, simulation_data)

  }

  simulation_data_df$type <- factor(simulation_data_df$type, levels = c('500 obs', '250 obs', '125 obs','60 obs'))
  #if(prior_sim == 'flat'){prior_sim = 'Flat'}
  #if(prior_sim == 'jeffreys'){prior_sim = 'Jeffreys'}
  #if(prior_sim == 'mdi'){prior_sim = 'MDI'}

  #Plot posterior density
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

  ggsave(here::here('05-week','figures_2','posterior_mmax_mean_different_obs', file_name),
         plot = plot_post,
         units = "px",
         width = 1197,
         height = 683,
         dpi = 300)

  #Plot 0.5 quantile
  quantile<- 0.5
  true_quantile <- quantile_posterior_calculation(mmax = mmax_sim ,mean = mean_sim, quantile =quantile, u = u_sim)
  plot_quantile <- ggplot(simulation_data_df, aes(x = quantile_0.5, colour = type, fill = type)) +
    geom_density(alpha = 0.4)+
    geom_segment(aes(x = true_quantile, y = 0, xend = true_quantile, yend = Inf, linetype = "True quantile"), color = 'black')+
    ggtitle(paste('Posterior distribution of quantile', quantile, sep = ' '))+
    theme_classic()+
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(family = 'sans'),
          axis.text.y = element_text(family = 'sans'),
          title = element_text(family = 'sans', size = 8),
          legend.title = element_blank(),
          legend.text = element_text(family = 'sans', size = 7),
          #panel.grid.major = element_blank(),
          #panel.grid.minor = element_blank(),
          axis.line.y = element_line(),
          axis.ticks.y = element_line(),
          axis.line.x = element_line(),
          axis.ticks.x = element_line())+
    #legend.direction='vertical')+
    #legend.position = 'right')+
    expand_limits(y = 0) +
    coord_cartesian(expand = FALSE, clip = "off")


  #save plot
  file_name <- paste('posterior_0.5_quantile',
                     prior_sim, 'prior',
                     mmax_sim, 'mmax',
                     mean_sim, 'mean.png', sep = '-')

  ggsave(here::here('05-week', 'figures_2','posterior_mmax_mean_different_obs', file_name),
         plot = plot_quantile,
         units = "px",
         width = 1197,
         height = 683,
         dpi = 300)

  #Plot 0.75 quantile
  quantile<- 0.75
  true_quantile <- quantile_posterior_calculation(mmax = mmax_sim ,mean = mean_sim, quantile =quantile, u = u_sim)
  plot_quantile <- ggplot(simulation_data_df, aes(x = quantile_0.75, colour = type, fill = type)) +
    geom_density(alpha = 0.4)+
    geom_segment(aes(x = true_quantile, y = 0, xend = true_quantile, yend = Inf, linetype = "True quantile"), color = 'black')+
    ggtitle(paste('Posterior distribution of quantile', quantile, sep = ' '))+
    theme_classic()+
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(family = 'sans'),
          axis.text.y = element_text(family = 'sans'),
          title = element_text(family = 'sans', size = 8),
          legend.title = element_blank(),
          legend.text = element_text(family = 'sans', size = 7),
          #panel.grid.major = element_blank(),
          #panel.grid.minor = element_blank(),
          axis.line.y = element_line(),
          axis.ticks.y = element_line(),
          axis.line.x = element_line(),
          axis.ticks.x = element_line())+
    #legend.direction='vertical')+
    #legend.position = 'right')+
    expand_limits(y = 0) +
    coord_cartesian(expand = FALSE, clip = "off")


  #save plot
  file_name <- paste('posterior_0.75_quantile',
                     prior_sim, 'prior',
                     mmax_sim, 'mmax',
                     mean_sim, 'mean.png', sep = '-')

  ggsave(here::here('05-week', 'figures_2','posterior_mmax_mean_different_obs', file_name),
         plot = plot_quantile,
         units = "px",
         width = 1197,
         height = 683,
         dpi = 300)


  #Plot 0.9 quantile
  quantile<- 0.9
  true_quantile <- quantile_posterior_calculation(mmax = mmax_sim ,mean = mean_sim, quantile =quantile, u = u_sim)
  plot_quantile <- ggplot(simulation_data_df, aes(x = quantile_0.9, colour = type, fill = type)) +
    geom_density(alpha = 0.4)+
    geom_segment(aes(x = true_quantile, y = 0, xend = true_quantile, yend = Inf, linetype = "True quantile"), color = 'black')+
    ggtitle(paste('Posterior distribution of quantile', quantile, sep = ' '))+
    theme_classic()+
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(family = 'sans'),
          axis.text.y = element_text(family = 'sans'),
          title = element_text(family = 'sans', size = 8),
          legend.title = element_blank(),
          legend.text = element_text(family = 'sans', size = 7),
          #panel.grid.major = element_blank(),
          #panel.grid.minor = element_blank(),
          axis.line.y = element_line(),
          axis.ticks.y = element_line(),
          axis.line.x = element_line(),
          axis.ticks.x = element_line())+
    #legend.direction='vertical')+
    #legend.position = 'right')+
    expand_limits(y = 0) +
    coord_cartesian(expand = FALSE, clip = "off")


  #save plot
  file_name <- paste('posterior_0.9_quantile',
                     prior_sim, 'prior',
                     mmax_sim, 'mmax',
                     mean_sim, 'mean.png', sep = '-')

  ggsave(here::here('05-week', 'figures_2','posterior_mmax_mean_different_obs', file_name),
         plot = plot_quantile,
         units = "px",
         width = 1197,
         height = 683,
         dpi = 300)

  #Plot 0.95 quantile
  quantile<- 0.95
  true_quantile <- quantile_posterior_calculation(mmax = mmax_sim ,mean = mean_sim, quantile =quantile, u = u_sim)
  plot_quantile <- ggplot(simulation_data_df, aes(x = quantile_0.95, colour = type, fill = type)) +
    geom_density(alpha = 0.4)+
    geom_segment(aes(x = true_quantile, y = 0, xend = true_quantile, yend = Inf, linetype = "True quantile"), color = 'black')+
    ggtitle(paste('Posterior distribution of quantile', quantile, sep = ' '))+
    theme_classic()+
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(family = 'sans'),
          axis.text.y = element_text(family = 'sans'),
          title = element_text(family = 'sans', size = 8),
          legend.title = element_blank(),
          legend.text = element_text(family = 'sans', size = 7),
          #panel.grid.major = element_blank(),
          #panel.grid.minor = element_blank(),
          axis.line.y = element_line(),
          axis.ticks.y = element_line(),
          axis.line.x = element_line(),
          axis.ticks.x = element_line())+
    #legend.direction='vertical')+
    #legend.position = 'right')+
    expand_limits(y = 0) +
    coord_cartesian(expand = FALSE, clip = "off")


  #save plot
  file_name <- paste('posterior_0.95_quantile',
                     prior_sim, 'prior',
                     mmax_sim, 'mmax',
                     mean_sim, 'mean.png', sep = '-')

  ggsave(here::here('05-week', 'figures_2','posterior_mmax_mean_different_obs', file_name),
         plot = plot_quantile,
         units = "px",
         width = 1197,
         height = 683,
         dpi = 300)
}


############################
#Analysis and visualization 3
############################
mmax_sim <- 7.6
mean_sim <- 2.1
u_sim <- 1.45
b_value_sim <- 1.8
epsilon_sim <- 1.5
upper_mmax_sim <- 8.5
alpha_sim_1 <- 1
beta_sim_1 <- 0.5
alpha_sim_2 <- 1
beta_sim_2 <- 0.5
prior_sim_vector <- c('unif-unif', 'unif-gamma', 'flat-gamma', 'gamma-gamma')
#n_data_sim_vector <- c(60, 125, 250, 500)
n_data_sim_vector <- c(25, 35, 45, 55)

for (j in 1:length(n_data_sim_vector)){
  n_data_sim <- n_data_sim_vector[j]
  #obtain data frame with results from different simulations
  simulation_data_df <- data.frame()
  for (i in 1:length(prior_sim_vector)){
    prior_sim <- prior_sim_vector[i]
    file_name_sim <- paste('mmax', mmax_sim,
                           'mean', mean_sim,
                           'u', u_sim,
                           prior_sim, 'prior',
                           n_data_sim, 'n_data', sep = '_')

    simulation_data <- read.csv(here::here('05-week','outputs_2',file_name_sim),
                                colClasses=c("NULL", NA, NA), col.names = c('','mmax', 'mean'))

    #add a column to indicate the number of observations used for the simulation
    simulation_data$type <- rep(paste(prior_sim, 'prior', sep = ' '), length(simulation_data$mean))

    #add a column to compute the 0.5, 0.75, 0.9 and 0.95 quantiles
    simulation_data$quantile_0.5 <- quantile_posterior_calculation(mmax = simulation_data$mmax, mean = simulation_data$mean, u = u_sim, quantile = 0.5)
    simulation_data$quantile_0.75 <- quantile_posterior_calculation(mmax = simulation_data$mmax, mean = simulation_data$mean, u = u_sim, quantile = 0.75)
    simulation_data$quantile_0.9 <- quantile_posterior_calculation(mmax = simulation_data$mmax, mean = simulation_data$mean, u = u_sim, quantile = 0.9)
    simulation_data$quantile_0.95 <- quantile_posterior_calculation(mmax = simulation_data$mmax, mean = simulation_data$mean, u = u_sim, quantile = 0.95)

    simulation_data_df <- rbind(simulation_data_df, simulation_data)

  }

  #simulation_data_df$type <- factor(simulation_data_df$type, levels = c('500 obs', '250 obs', '125 obs','60 obs'))
  #if(prior_sim == 'flat'){prior_sim = 'Flat'}
  #if(prior_sim == 'jeffreys'){prior_sim = 'Jeffreys'}
  #if(prior_sim == 'mdi'){prior_sim = 'MDI'}

  #Plot posterior density
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

  ggsave(here::here('05-week','figures_2','posterior_mmax_mean_different_obs', file_name),
         plot = plot_post,
         units = "px",
         width = 1197,
         height = 683,
         dpi = 300)

  #Plot 0.5 quantile
  quantile<- 0.5
  true_quantile <- quantile_posterior_calculation(mmax = mmax_sim ,mean = mean_sim, quantile =quantile, u = u_sim)
  plot_quantile <- ggplot(simulation_data_df, aes(x = quantile_0.5, colour = type, fill = type)) +
    geom_density(alpha = 0.4)+
    geom_segment(aes(x = true_quantile, y = 0, xend = true_quantile, yend = Inf, linetype = "True quantile"), color = 'black')+
    ggtitle(paste('Posterior distribution of quantile', quantile, sep = ' '))+
    theme_classic()+
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(family = 'sans'),
          axis.text.y = element_text(family = 'sans'),
          title = element_text(family = 'sans', size = 8),
          legend.title = element_blank(),
          legend.text = element_text(family = 'sans', size = 7),
          #panel.grid.major = element_blank(),
          #panel.grid.minor = element_blank(),
          axis.line.y = element_line(),
          axis.ticks.y = element_line(),
          axis.line.x = element_line(),
          axis.ticks.x = element_line())+
    #legend.direction='vertical')+
    #legend.position = 'right')+
    expand_limits(y = 0) +
    coord_cartesian(expand = FALSE, clip = "off")


  #save plot
  file_name <- paste('posterior_0.5_quantile',
                     prior_sim, 'prior',
                     mmax_sim, 'mmax',
                     mean_sim, 'mean.png', sep = '-')

  ggsave(here::here('05-week', 'figures_2','posterior_mmax_mean_different_obs', file_name),
         plot = plot_quantile,
         units = "px",
         width = 1197,
         height = 683,
         dpi = 300)

  #Plot 0.75 quantile
  quantile<- 0.75
  true_quantile <- quantile_posterior_calculation(mmax = mmax_sim ,mean = mean_sim, quantile =quantile, u = u_sim)
  plot_quantile <- ggplot(simulation_data_df, aes(x = quantile_0.75, colour = type, fill = type)) +
    geom_density(alpha = 0.4)+
    geom_segment(aes(x = true_quantile, y = 0, xend = true_quantile, yend = Inf, linetype = "True quantile"), color = 'black')+
    ggtitle(paste('Posterior distribution of quantile', quantile, sep = ' '))+
    theme_classic()+
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(family = 'sans'),
          axis.text.y = element_text(family = 'sans'),
          title = element_text(family = 'sans', size = 8),
          legend.title = element_blank(),
          legend.text = element_text(family = 'sans', size = 7),
          #panel.grid.major = element_blank(),
          #panel.grid.minor = element_blank(),
          axis.line.y = element_line(),
          axis.ticks.y = element_line(),
          axis.line.x = element_line(),
          axis.ticks.x = element_line())+
    #legend.direction='vertical')+
    #legend.position = 'right')+
    expand_limits(y = 0) +
    coord_cartesian(expand = FALSE, clip = "off")


  #save plot
  file_name <- paste('posterior_0.75_quantile',
                     prior_sim, 'prior',
                     mmax_sim, 'mmax',
                     mean_sim, 'mean.png', sep = '-')

  ggsave(here::here('05-week', 'figures_2','posterior_mmax_mean_different_obs', file_name),
         plot = plot_quantile,
         units = "px",
         width = 1197,
         height = 683,
         dpi = 300)


  #Plot 0.9 quantile
  quantile<- 0.9
  true_quantile <- quantile_posterior_calculation(mmax = mmax_sim ,mean = mean_sim, quantile =quantile, u = u_sim)
  plot_quantile <- ggplot(simulation_data_df, aes(x = quantile_0.9, colour = type, fill = type)) +
    geom_density(alpha = 0.4)+
    geom_segment(aes(x = true_quantile, y = 0, xend = true_quantile, yend = Inf, linetype = "True quantile"), color = 'black')+
    ggtitle(paste('Posterior distribution of quantile', quantile, sep = ' '))+
    theme_classic()+
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(family = 'sans'),
          axis.text.y = element_text(family = 'sans'),
          title = element_text(family = 'sans', size = 8),
          legend.title = element_blank(),
          legend.text = element_text(family = 'sans', size = 7),
          #panel.grid.major = element_blank(),
          #panel.grid.minor = element_blank(),
          axis.line.y = element_line(),
          axis.ticks.y = element_line(),
          axis.line.x = element_line(),
          axis.ticks.x = element_line())+
    #legend.direction='vertical')+
    #legend.position = 'right')+
    expand_limits(y = 0) +
    coord_cartesian(expand = FALSE, clip = "off")


  #save plot
  file_name <- paste('posterior_0.9_quantile',
                     prior_sim, 'prior',
                     mmax_sim, 'mmax',
                     mean_sim, 'mean.png', sep = '-')

  ggsave(here::here('05-week', 'figures_2','posterior_mmax_mean_different_obs', file_name),
         plot = plot_quantile,
         units = "px",
         width = 1197,
         height = 683,
         dpi = 300)

  #Plot 0.95 quantile
  quantile<- 0.95
  true_quantile <- quantile_posterior_calculation(mmax = mmax_sim ,mean = mean_sim, quantile =quantile, u = u_sim)
  plot_quantile <- ggplot(simulation_data_df, aes(x = quantile_0.95, colour = type, fill = type)) +
    geom_density(alpha = 0.4)+
    geom_segment(aes(x = true_quantile, y = 0, xend = true_quantile, yend = Inf, linetype = "True quantile"), color = 'black')+
    ggtitle(paste('Posterior distribution of quantile', quantile, sep = ' '))+
    theme_classic()+
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(family = 'sans'),
          axis.text.y = element_text(family = 'sans'),
          title = element_text(family = 'sans', size = 8),
          legend.title = element_blank(),
          legend.text = element_text(family = 'sans', size = 7),
          #panel.grid.major = element_blank(),
          #panel.grid.minor = element_blank(),
          axis.line.y = element_line(),
          axis.ticks.y = element_line(),
          axis.line.x = element_line(),
          axis.ticks.x = element_line())+
    #legend.direction='vertical')+
    #legend.position = 'right')+
    expand_limits(y = 0) +
    coord_cartesian(expand = FALSE, clip = "off")


  #save plot
  file_name <- paste('posterior_0.95_quantile',
                     prior_sim, 'prior',
                     mmax_sim, 'mmax',
                     mean_sim, 'mean.png', sep = '-')

  ggsave(here::here('05-week', 'figures_2','posterior_mmax_mean_different_obs', file_name),
         plot = plot_quantile,
         units = "px",
         width = 1197,
         height = 683,
         dpi = 300)
}
