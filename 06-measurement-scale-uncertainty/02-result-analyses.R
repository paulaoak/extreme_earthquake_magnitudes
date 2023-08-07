###################################
#READ DATA AND COMPUTE QUANTILES
###################################


library(ggplot2)
#library(urbnthemes)


#function to compute quantile from xi and sigma
quantile_posterior_calculation_penultimate <- function(xi, sigma, quantile, u){
  p <- 1 - quantile

  quantiles_posterior <- u + sigma/xi * (p^(-xi) - 1)
  return(quantiles_posterior)
}

#function to compute quantile from xi, sigma and scale
quantile_posterior_calculation_no_penultimate <- function(xi, sigma, lamdba, quantile, u){
  p <- 1 - quantile

  u_mod <- (u^lamdba-1)/lamdba
  quantiles_posterior_transform <- u_mod + sigma/xi * (p^(-xi) - 1)
  quantiles_posterior <- (quantiles_posterior_transform * lamdba + 1)^(1/lamdba)

  return(quantiles_posterior)
}


#function to compute quantile from xi and sigma
quantile_penultimate <- function(xi, sigma, quantile, u, transform){
  p <- 1 - quantile
  u <- ifelse(transform, log(u), u)
  quantiles_posterior <- u + sigma/xi * (p^(-xi) - 1)
  return(quantiles_posterior)
}

#function to compute quantile from xi, sigma and scale
quantile_no_penultimate <- function(xi, sigma, lamdba, quantile, u, transform){
  p <- 1 - quantile
  u <- ifelse(transform, log(u), u)
  u_mod <- (u^lamdba-1)/lamdba
  quantiles_posterior_transform <- u_mod + sigma/xi * (p^(-xi) - 1)
  quantiles_posterior <- (quantiles_posterior_transform * lamdba + 1)^(1/lamdba)

  return(quantiles_posterior)
}

############################
#Analysis and visualization 1
############################
xi_sim <- -0.084
sigma_sim <- 0.48
u_sim <- 1.45
#n_data_sim_vector <- c(60, 75, 125, 250, 500)
n_data_sim <- 75
#prior_sim <- 'flat'
prior_sim <- 'mdi'
scale_transformation_sim_vector <- c(TRUE, FALSE)

#obtain data frame with results from different simulations
simulation_data_df <- data.frame()

#READ FILES

#######################
#File 1
#######################
file_name_sim_penultimate_true <- paste('new_simulation_penultimate', 'xi', xi_sim,
                                     'sigma', sigma_sim, 'u', u_sim,
                                     prior_sim, 'prior', n_data_sim, 'n_data',
                                     'transform_data',scale_transformation_sim_vector[1], sep = '_')

simulation_data <- read.csv(here::here('06-measurement-scale-uncertainty','outputs', file_name_sim_penultimate_true),
                            colClasses=c("NULL", NA, NA, NA), col.names = c('','sigma', 'xi', 'lambda'))

simulation_data$transformation <- rep(TRUE, length(simulation_data$sigma))
simulation_data$type <- rep('penultimate', length(simulation_data$sigma))
simulation_data$class <- rep('Penultimate, Log data', length(simulation_data$sigma))
simulation_data_df <- rbind(simulation_data_df, simulation_data)

u_vector <- rep(log(u_sim), length(simulation_data$sigma))
#######################
#File 2
#######################
file_name_sim_2_true <- paste('new_simulation', 'xi', xi_sim,
                        'sigma', sigma_sim, 'u', u_sim,
                        prior_sim, 'prior', n_data_sim, 'n_data',
                        'transform_data',scale_transformation_sim_vector[1], sep = '_')

simulation_data <- read.csv(here::here('06-measurement-scale-uncertainty','outputs', file_name_sim_2_true),
                            colClasses=c("NULL", NA, NA, NA), col.names = c('','sigma', 'xi', 'lambda'))

simulation_data$transformation <- rep(TRUE, length(simulation_data$sigma))
simulation_data$type <- rep('no_penultimate', length(simulation_data$sigma))
simulation_data$class <- rep('No penultimate, Log data', length(simulation_data$sigma))
simulation_data_df <- rbind(simulation_data_df, simulation_data)

u_vector <- c(u_vector, rep(log(u_sim), length(simulation_data$sigma)))
#######################
#File 3
#######################
file_name_sim_penultimate_false <- paste('new_simulation_penultimate', 'xi', xi_sim,
                                        'sigma', sigma_sim, 'u', u_sim,
                                        prior_sim, 'prior', n_data_sim, 'n_data',
                                        'transform_data',scale_transformation_sim_vector[2], sep = '_')

simulation_data <- read.csv(here::here('06-measurement-scale-uncertainty','outputs', file_name_sim_penultimate_false),
                            colClasses=c("NULL", NA, NA, NA), col.names = c('','sigma', 'xi', 'lambda'))

simulation_data$transformation <- rep(FALSE, length(simulation_data$sigma))
simulation_data$type <- rep('penultimate', length(simulation_data$sigma))
simulation_data$class <- rep('Penultimate', length(simulation_data$sigma))
simulation_data_df <- rbind(simulation_data_df, simulation_data)

u_vector <- c(u_vector, rep(u_sim, length(simulation_data$sigma)))
#######################
#File 4
#######################
file_name_sim_2_false <- paste('new_simulation', 'xi', xi_sim,
                              'sigma', sigma_sim, 'u', u_sim,
                              prior_sim, 'prior', n_data_sim, 'n_data',
                              'transform_data',scale_transformation_sim_vector[2], sep = '_')

simulation_data <- read.csv(here::here('06-measurement-scale-uncertainty','outputs', file_name_sim_2_false),
                            colClasses=c("NULL", NA, NA, NA), col.names = c('','sigma', 'xi', 'lambda'))

simulation_data$transformation <- rep(FALSE, length(simulation_data$sigma))
simulation_data$type <- rep('no_penultimate', length(simulation_data$sigma))
simulation_data$class <- rep('No penultimate', length(simulation_data$sigma))
simulation_data_df <- rbind(simulation_data_df, simulation_data)

u_vector <- c(u_vector, rep(u_sim, length(simulation_data$sigma)))

#######################
#File 5
#######################
file_name_sim_no_false <- paste('new_simulation_NO', 'xi', xi_sim,
                                'sigma', sigma_sim, 'u', u_sim,
                                prior_sim, 'prior', n_data_sim, 'n_data',
                                'transform_data',scale_transformation_sim_vector[2], sep = '_')

simulation_data <- read.csv(here::here('06-measurement-scale-uncertainty','outputs', file_name_sim_no_false),
                            colClasses=c("NULL", NA, NA, NA), col.names = c('','sigma', 'xi', 'lambda'))

simulation_data$transformation <- rep(FALSE, length(simulation_data$sigma))
simulation_data$type <- rep('no_scale', length(simulation_data$sigma))
simulation_data$class <- rep('No scale', length(simulation_data$sigma))
simulation_data_df <- rbind(simulation_data_df, simulation_data)
#u_vector <- c()
u_vector <- c(u_vector, rep(u_sim, length(simulation_data$sigma)))
#######################
#File 6
#######################
file_name_sim_no_true <- paste('new_simulation_NO', 'xi', xi_sim,
                               'sigma', sigma_sim, 'u', u_sim,
                               prior_sim, 'prior', n_data_sim, 'n_data',
                               'transform_data',scale_transformation_sim_vector[1], sep = '_')

simulation_data <- read.csv(here::here('06-measurement-scale-uncertainty','outputs', file_name_sim_no_true),
                            colClasses=c("NULL", NA, NA, NA), col.names = c('','sigma', 'xi', 'lambda'))

simulation_data$transformation <- rep(TRUE, length(simulation_data$sigma))
simulation_data$type <- rep('no_scale', length(simulation_data$sigma))
simulation_data$class <- rep('No scale, Log data', length(simulation_data$sigma))
simulation_data_df <- rbind(simulation_data_df, simulation_data)

u_vector <- c(u_vector, rep(log(u_sim), length(simulation_data$sigma)))
#u_vector <- c(rep(log(u_sim), length(simulation_data$sigma)))
#################################
# Compute posterior quantiles
#################################
#add a column to compute the 0.5, 0.75, 0.9 and 0.95 quantiles
#simulation_data_df$quantile_0.5 <- ifelse(simulation_data_df$type=='penultimate',
#                                          quantile_posterior_calculation_penultimate(sigma = simulation_data_df$sigma, xi = simulation_data_df$xi, u = u_vector, quantile = 0.5),
#                                          quantile_posterior_calculation_no_penultimate(sigma = simulation_data_df$sigma, xi = simulation_data_df$xi, u = u_vector, quantile = 0.5, lamdba = simulation_data_df$lambda))

simulation_data_df$quantile_0.5 <- quantile_posterior_calculation_penultimate(sigma = simulation_data_df$sigma, xi = simulation_data_df$xi, u = u_vector, quantile = 0.5)
simulation_data_df$quantile_0.5 <- ifelse(simulation_data_df$transformation, exp(simulation_data_df$quantile_0.5), simulation_data_df$quantile_0.5)

#simulation_data_df$quantile_0.75 <- ifelse(simulation_data_df$type=='penultimate',
#                                           quantile_posterior_calculation_penultimate(sigma = simulation_data_df$sigma, xi = simulation_data_df$xi, u = u_vector, quantile = 0.75),
#                                           quantile_posterior_calculation_no_penultimate(sigma = simulation_data_df$sigma, xi = simulation_data_df$xi, u = u_vector, quantile = 0.75, lamdba = simulation_data_df$lambda))

simulation_data_df$quantile_0.75 <- quantile_posterior_calculation_penultimate(sigma = simulation_data_df$sigma, xi = simulation_data_df$xi, u = u_vector, quantile = 0.75)
simulation_data_df$quantile_0.75 <- ifelse(simulation_data_df$transformation, exp(simulation_data_df$quantile_0.75), simulation_data_df$quantile_0.75)

simulation_data_df$quantile_0.9 <- ifelse(simulation_data_df$type == 'no_penultimate',
                                          quantile_no_penultimate(sigma = simulation_data_df$sigma, xi = simulation_data_df$xi, u = u_sim, quantile = 0.9, lamdba = simulation_data_df$lambda, transform = simulation_data_df$transformation),
                                          quantile_penultimate(sigma = simulation_data_df$sigma, xi = simulation_data_df$xi, u = u_sim, quantile = 0.9, transform = simulation_data_df$transformation))

#simulation_data_df$quantile_0.9 <- quantile_posterior_calculation_penultimate(sigma = simulation_data_df$sigma, xi = simulation_data_df$xi, u = u_vector, quantile = 0.9)
simulation_data_df$quantile_0.9 <- ifelse(simulation_data_df$transformation, exp(simulation_data_df$quantile_0.9), simulation_data_df$quantile_0.9)

#simulation_data_df$quantile_0.95 <- ifelse(simulation_data_df$type == 'penultimate',
#                                           quantile_posterior_calculation_penultimate(sigma = simulation_data_df$sigma, xi = simulation_data_df$xi, u = u_vector, quantile = 0.95),
#                                           quantile_posterior_calculation_no_penultimate(sigma = simulation_data_df$sigma, xi = simulation_data_df$xi, u = u_vector, quantile = 0.95, lamdba = simulation_data_df$lambda))

simulation_data_df$quantile_0.95 <- quantile_posterior_calculation_penultimate(sigma = simulation_data_df$sigma, xi = simulation_data_df$xi, u = u_vector, quantile = 0.95)
simulation_data_df$quantile_0.95 <- ifelse(simulation_data_df$transformation, exp(simulation_data_df$quantile_0.95), simulation_data_df$quantile_0.95)

#simulation_data_df$type <- factor(simulation_data_df$type, levels = c('500 obs', '250 obs', '125 obs','60 obs'))
#if(prior_sim == 'flat'){prior_sim = 'Flat'}
#if(prior_sim == 'jeffreys'){prior_sim = 'Jeffreys'}
#if(prior_sim == 'mdi'){prior_sim = 'MDI'}


#Plot 0.5 quantile
quantile<- 0.5
true_quantile <- quantile_posterior_calculation_penultimate(sigma = sigma_sim ,xi = xi_sim, quantile =quantile, u = u_sim)
plot_quantile <- ggplot(simulation_data_df, aes(x = quantile_0.5, colour = class, fill = class)) +
  geom_density(alpha = 0.25)+
  geom_segment(aes(x = true_quantile, y = 0, xend = true_quantile, yend = Inf, linetype = "True quantile"), color = 'black')+
  #ggtitle(paste('Posterior distribution of quantile', quantile, sep = ' '))+
  xlab(paste(quantile, 'quantile zoom', sep = ' '))+
  theme_classic()+
  theme(#axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(family = 'sans'),
        axis.text.y = element_text(family = 'sans'),
        title = element_text(family = 'sans', size = 8),
        legend.title = element_blank(),
        legend.text = element_text(family = 'sans', size = 8),
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        axis.line.y = element_line(),
        axis.ticks.y = element_line(),
        axis.line.x = element_line(),
        axis.ticks.x = element_line())+
  #legend.direction='vertical')+
  #legend.position = 'right')+
  expand_limits(y = 0) +
  coord_cartesian(expand = FALSE, clip = "off")+
  xlim(1.6,2.2)


#save plot
file_name <- paste('posterior_0.5_quantile_zoom',
                   prior_sim, 'prior',
                   sigma_sim, 'sigma',
                   xi_sim, 'xi.png', sep = '-')

ggsave(here::here('06-measurement-scale-uncertainty', 'figures', file_name),
       plot = plot_quantile,
       units = "px",
       width = 1197,
       height = 683,
       dpi = 300)

#Plot 0.75 quantile
quantile<- 0.75
true_quantile <- quantile_posterior_calculation_penultimate(sigma = sigma_sim ,xi = xi_sim, quantile =quantile, u = u_sim)
plot_quantile <- ggplot(simulation_data_df, aes(x = quantile_0.75, colour = class, fill = class)) +
  geom_density(alpha = 0.25)+
  geom_segment(aes(x = true_quantile, y = 0, xend = true_quantile, yend = Inf, linetype = "True quantile"), color = 'black')+
  #ggtitle(paste('Posterior distribution of quantile', quantile, sep = ' '))+
  xlab(paste(quantile, 'quantile', sep = ' '))+
  theme_classic()+
  theme(#axis.title.x = element_blank(),
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
  coord_cartesian(expand = FALSE, clip = "off")+
  xlim(1.4, 4)


#save plot
file_name <- paste('posterior_0.75_quantile',
                   prior_sim, 'prior',
                   sigma_sim, 'sigma',
                   xi_sim, 'xi.png', sep = '-')

ggsave(here::here('06-measurement-scale-uncertainty', 'figures', file_name),
       plot = plot_quantile,
       units = "px",
       width = 1197,
       height = 683,
       dpi = 300)


#Plot 0.9 quantile
quantile<- 0.9
true_quantile <- quantile_posterior_calculation_penultimate(sigma = sigma_sim ,xi = xi_sim, quantile =quantile, u = u_sim)
plot_quantile <- ggplot(simulation_data_df, aes(x = quantile_0.9, colour = class, fill = class)) +
  geom_density(alpha = 0.25)+
  geom_segment(aes(x = true_quantile, y = 0, xend = true_quantile, yend = Inf, linetype = "True quantile"), color = 'black')+
  #ggtitle(paste('Posterior distribution of quantile', quantile, sep = ' '))+
  xlab(paste(quantile, 'quantile', sep = ' '))+
  theme_classic()+
  theme(#axis.title.x = element_blank(),
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
  coord_cartesian(expand = FALSE, clip = "off")+
  xlim(1.4, 5)


#save plot
file_name <- paste('posterior_0.9_quantile',
                   prior_sim, 'prior',
                   sigma_sim, 'sigma',
                   xi_sim, 'xi.png', sep = '-')

ggsave(here::here('06-measurement-scale-uncertainty', 'figures', file_name),
       plot = plot_quantile,
       units = "px",
       width = 1197,
       height = 683,
       dpi = 300)

#Plot 0.95 quantile
quantile<- 0.95
true_quantile <- quantile_posterior_calculation_penultimate(sigma = sigma_sim ,xi = xi_sim, quantile =quantile, u = u_sim)
plot_quantile <- ggplot(simulation_data_df, aes(x = quantile_0.95, colour = class, fill = class)) +
  geom_density(alpha = 0.25)+
  geom_segment(aes(x = true_quantile, y = 0, xend = true_quantile, yend = Inf, linetype = "True quantile"), color = 'black')+
  #ggtitle(paste('Posterior distribution of quantile', quantile, sep = ' '))+
  xlab(paste(quantile, 'quantile', sep = ' ')) +
  theme_classic()+
  theme(#axis.title.x = element_blank(),
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
  coord_cartesian(expand = FALSE, clip = "off")+
  xlim(1.25, 5.5)


#save plot
file_name <- paste('posterior_0.95_quantile',
                   prior_sim, 'prior',
                   sigma_sim, 'sigma',
                   xi_sim, 'xi.png', sep = '-')

ggsave(here::here('06-measurement-scale-uncertainty', 'figures', file_name),
       plot = plot_quantile,
       units = "px",
       width = 1197,
       height = 683,
       dpi = 300)
