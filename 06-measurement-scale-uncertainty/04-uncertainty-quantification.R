###################################
#READ DATA AND COMPUTE QUANTILES
###################################


library(ggplot2)

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

############################
#Analysis and visualization 1
############################
xi_sim <- -0.084
sigma_sim <- 0.48
u_sim <- 1.45
#n_data_sim_vector <- c(60, 75, 125, 250, 500)
n_data_sim <- 75
prior_sim <- 'flat'
#prior_sim <- 'mdi'
scale_transformation_sim_vector <- c(TRUE, FALSE)
quantile_value <- 0.95

#obtain data frame with results from different simulations
simulation_data_df <- data.frame()

#READ FILES

#######################
#File 1
#######################
file_name_sim_penultimate_true <- paste('simulation_penultimate_linear', 'xi', xi_sim,
                                        'sigma', sigma_sim, 'u', u_sim,
                                        prior_sim, 'prior', n_data_sim, 'n_data',
                                        'transform_data',scale_transformation_sim_vector[1], sep = '_')

simulation_data <- read.csv(here::here('06-measurement-scale-uncertainty','outputs', file_name_sim_penultimate_true),
                            colClasses=c("NULL", NA, NA, NA), col.names = c('','sigma', 'xi', 'lambda'))

simulation_data$transformation <- rep(TRUE, length(simulation_data$sigma))
u_vector <- rep(log(u_sim), length(simulation_data$sigma))

simulation_data$quantile_0.75 <- quantile_posterior_calculation_penultimate(sigma = simulation_data$sigma, xi = simulation_data$xi, u = u_vector, quantile = quantile_value)
simulation_data$quantile_0.75 <- ifelse(simulation_data$transformation, exp(simulation_data$quantile_0.75), simulation_data$quantile_0.75)

df_quantile_75_penultimate_true <- approxfun(density(simulation_data$quantile_0.75))

plot(density(simulation_data$quantile_0.75))
x_plot<-seq(0, 7, by = 0.001)
lines(x_plot ,df_quantile_75_penultimate_true(x_plot), type = 'l', col='red')
true_penultimate <- df_quantile_75_penultimate_true(x_plot)
true_penultimate <- ifelse(is.na(true_penultimate), 0, true_penultimate)
simulation_data_aux<-data.frame(x_plot = x_plot, y_plot = true_penultimate,
                                    type = rep('Single model, Log data', length(x_plot)))

simulation_data_df <- rbind(simulation_data_df, simulation_data_aux)
#######################
#File 2
#######################
c_hat_sig_c <- simulation_mcmc_scale_penultimate_uq_c_estimates(xi = xi_sim, sigma = sigma_sim,
                                                 xi_init_min = -0.5, sigma_init_min = 0,
                                                 xi_init_max = 1, sigma_init_max = 4,
                                                 n_starting_points = 5, u = u_sim,
                                                 n_data = n_data_sim, n_iter = 1e4, n_burn = 1e3,
                                                 prior = prior_sim,
                                                 sd_xi = c(0.2, 0.15, 0.1, 0.05), sd_sigma = c(0.2, 0.15, 0.1, 0.05),
                                                 min_lambda = 0, max_lambda = 2.1, step_lambda = 0.1,
                                                 min_xi = -0.35, max_xi = 0.35, step_xi = 0.05,
                                                 scale_transformation = scale_transformation_sim_vector[1])

c_hat <- c_hat_sig_c[1]
sig_c <- c_hat_sig_c[2]
file_name_sim_2_true <- paste('UQ_simulation_penultimate', 'xi', xi_sim,
                              'sigma', sigma_sim, 'u', u_sim,
                              prior_sim, 'prior', n_data_sim, 'n_data',
                              'transform_data',scale_transformation_sim_vector[1], sep = '_')

simulation_data <- read.csv(here::here('06-measurement-scale-uncertainty','outputs', file_name_sim_2_true),
                            colClasses=c("NULL", NA, NA, NA, NA), col.names = c('','sigma', 'xi', 'lambda', 'c_value'))

simulation_data$transformation <- rep(TRUE, length(simulation_data$sigma))
c_vector <- unique(simulation_data$c_value)
u_vector <- rep(log(u_sim), length(simulation_data$sigma))
simulation_data$quantile_0.75 <- quantile_posterior_calculation_penultimate(sigma = simulation_data$sigma, xi = simulation_data$xi, u = u_vector, quantile = quantile_value)
simulation_data$quantile_0.75 <- ifelse(simulation_data$transformation, exp(simulation_data$quantile_0.75), simulation_data$quantile_0.75)

df_quantile_75_c1 <- approxfun(density(filter(simulation_data, c_value == c_vector[1])$quantile_0.75))
df_quantile_75_c2 <- approxfun(density(filter(simulation_data, c_value == c_vector[2])$quantile_0.75))
df_quantile_75_c3 <- approxfun(density(filter(simulation_data, c_value == c_vector[3])$quantile_0.75))
df_quantile_75_c4 <- approxfun(density(filter(simulation_data, c_value == c_vector[4])$quantile_0.75))
df_quantile_75_c5 <- approxfun(density(filter(simulation_data, c_value == c_vector[5])$quantile_0.75))
df_quantile_75_c6 <- approxfun(density(filter(simulation_data, c_value == c_vector[6])$quantile_0.75))
df_quantile_75_c7 <- approxfun(density(filter(simulation_data, c_value == c_vector[7])$quantile_0.75))
#df_quantile_75_c8 <- approxfun(density(filter(simulation_data, c_value == c_vector[8])$quantile_0.75))

mixed_df_quantile_75_uq_true <- function(x){
  weights <- dnorm(c_vector, mean = c_hat, sd = sig_c)
  weights <- weights /sum(weights)
  df_quantile_75_c1(x) * weights[1] + df_quantile_75_c2(x) * weights[2] + df_quantile_75_c3(x) * weights[3] + df_quantile_75_c4(x) * weights[4] +
    df_quantile_75_c5(x) * weights[5] + df_quantile_75_c6(x) * weights[6] + df_quantile_75_c7(x) * weights[7] #+ df_quantile_75_c8(x) * weights[8]
}

plot(x_plot , mixed_df_quantile_75_uq_true(x_plot), type = 'l', col='red')
lines(x_plot ,df_quantile_75_penultimate_true(x_plot), type = 'l', col ='blue')
#lines(x_plot , df_quantile_75_c1(x_plot), type = 'l', col='black')
#lines(x_plot , df_quantile_75_c2(x_plot), type = 'l', col='black')
#lines(x_plot , df_quantile_75_c3(x_plot), type = 'l', col='black')
#lines(x_plot , df_quantile_75_c4(x_plot), type = 'l', col='black')
#lines(x_plot , df_quantile_75_c5(x_plot), type = 'l', col='black')
#lines(x_plot , df_quantile_75_c6(x_plot), type = 'l', col='black')
#lines(x_plot , df_quantile_75_c7(x_plot), type = 'l', col='black')
#lines(x_plot , df_quantile_75_c8(x_plot), type = 'l', col='black')
#simulation_data_df<-simulation_data

#u_vector <- c(u_vector, rep(log(u_sim), length(simulation_data$sigma)))
x_plot<-seq(0, 7, by = 0.001)
true_penultimate_uq <- mixed_df_quantile_75_uq_true(x_plot)
true_penultimate_uq <- ifelse(is.na(true_penultimate_uq), 0, true_penultimate_uq)
simulation_data_aux<-data.frame(x_plot = x_plot, y_plot = true_penultimate_uq,
                                     type = rep('Model avg., Log data', length(x_plot)))

simulation_data_df <- rbind(simulation_data_df, simulation_data_aux)
#######################
#File 3
#######################
file_name_sim_penultimate_false <- paste('simulation_penultimate_linear', 'xi', xi_sim,
                                         'sigma', sigma_sim, 'u', u_sim,
                                         prior_sim, 'prior', n_data_sim, 'n_data',
                                         'transform_data',scale_transformation_sim_vector[2], sep = '_')

simulation_data <- read.csv(here::here('06-measurement-scale-uncertainty','outputs', file_name_sim_penultimate_false),
                            colClasses=c("NULL", NA, NA, NA), col.names = c('','sigma', 'xi', 'lambda'))

simulation_data$transformation <- rep(FALSE, length(simulation_data$sigma))
u_vector <- rep(u_sim, length(simulation_data$sigma))

simulation_data$quantile_0.75 <- quantile_posterior_calculation_penultimate(sigma = simulation_data$sigma, xi = simulation_data$xi, u = u_vector, quantile = quantile_value)
simulation_data$quantile_0.75 <- ifelse(simulation_data$transformation, exp(simulation_data$quantile_0.75), simulation_data$quantile_0.75)

df_quantile_75_penultimate_false <- approxfun(density(simulation_data$quantile_0.75))

x_plot<-seq(0, 7, by = 0.001)
false_penultimate <- df_quantile_75_penultimate_false(x_plot)
false_penultimate <- ifelse(is.na(false_penultimate), 0, false_penultimate)
simulation_data_aux<-data.frame(x_plot = x_plot, y_plot = false_penultimate,
                                     type = rep('Single model', length(x_plot)))

simulation_data_df <- rbind(simulation_data_df, simulation_data_aux)
#######################
#File 4
#######################
c_hat_sig_c <- simulation_mcmc_scale_penultimate_uq_c_estimates(xi = xi_sim, sigma = sigma_sim,
                                                                xi_init_min = -0.5, sigma_init_min = 0,
                                                                xi_init_max = 1, sigma_init_max = 4,
                                                                n_starting_points = 5, u = u_sim,
                                                                n_data = n_data_sim, n_iter = 1e4, n_burn = 1e3,
                                                                prior = prior_sim,
                                                                sd_xi = c(0.2, 0.15, 0.1, 0.05), sd_sigma = c(0.2, 0.15, 0.1, 0.05),
                                                                min_lambda = 0, max_lambda = 2.1, step_lambda = 0.1,
                                                                min_xi = -0.35, max_xi = 0.35, step_xi = 0.05,
                                                                scale_transformation = scale_transformation_sim_vector[2])

c_hat <- c_hat_sig_c[1]
sig_c <- c_hat_sig_c[2]

file_name_sim_2_false <- paste('UQ_simulation_penultimate', 'xi', xi_sim,
                               'sigma', sigma_sim, 'u', u_sim,
                               prior_sim, 'prior', n_data_sim, 'n_data',
                               'transform_data',scale_transformation_sim_vector[2], sep = '_')

simulation_data <- read.csv(here::here('06-measurement-scale-uncertainty','outputs', file_name_sim_2_false),
                            colClasses=c("NULL", NA, NA, NA, NA), col.names = c('','sigma', 'xi', 'lambda', 'c_value'))

simulation_data$transformation <- rep(FALSE, length(simulation_data$sigma))

c_vector <- unique(simulation_data$c_value)
u_vector <- rep(u_sim, length(simulation_data$sigma))
simulation_data$quantile_0.75 <- quantile_posterior_calculation_penultimate(sigma = simulation_data$sigma, xi = simulation_data$xi, u = u_vector, quantile = quantile_value)
simulation_data$quantile_0.75 <- ifelse(simulation_data$transformation, exp(simulation_data$quantile_0.75), simulation_data$quantile_0.75)

df_quantile_75_c1 <- approxfun(density(filter(simulation_data, c_value == c_vector[1])$quantile_0.75))
df_quantile_75_c2 <- approxfun(density(filter(simulation_data, c_value == c_vector[2])$quantile_0.75))
df_quantile_75_c3 <- approxfun(density(filter(simulation_data, c_value == c_vector[3])$quantile_0.75))
df_quantile_75_c4 <- approxfun(density(filter(simulation_data, c_value == c_vector[4])$quantile_0.75))
df_quantile_75_c5 <- approxfun(density(filter(simulation_data, c_value == c_vector[5])$quantile_0.75))
df_quantile_75_c6 <- approxfun(density(filter(simulation_data, c_value == c_vector[6])$quantile_0.75))
df_quantile_75_c7 <- approxfun(density(filter(simulation_data, c_value == c_vector[7])$quantile_0.75))

mixed_df_quantile_75_uq_false <- function(x){
  weights <- dnorm(c_vector, mean = c_hat, sd = sig_c)
  weights <- weights /sum(weights)
  df_quantile_75_c1(x) * weights[1] + df_quantile_75_c2(x) * weights[2] + df_quantile_75_c3(x) * weights[3] + df_quantile_75_c4(x) * weights[4] +
    df_quantile_75_c5(x) * weights[5] + df_quantile_75_c6(x) * weights[6] + df_quantile_75_c7(x) * weights[7] #+ df_quantile_75_c8(x) * weights[8]
}

plot(x_plot , mixed_df_quantile_75_uq_false(x_plot), type = 'l', col='red')
lines(x_plot ,df_quantile_75_penultimate_false(x_plot), type = 'l', col ='blue')

x_plot<-seq(0, 7, by = 0.001)
false_penultimate_uq <- mixed_df_quantile_75_uq_false(x_plot)
false_penultimate_uq <- ifelse(is.na(false_penultimate_uq), 0, false_penultimate_uq)
simulation_data_aux<-data.frame(x_plot = x_plot, y_plot = false_penultimate_uq,
                                     type = rep('Model avg.', length(x_plot)))

simulation_data_df <- rbind(simulation_data_df, simulation_data_aux)

#u_vector <- c(u_vector, rep(u_sim, length(simulation_data$sigma)))


#Plot quantile
quantile<- quantile_value
true_quantile <- quantile_posterior_calculation_penultimate(sigma = sigma_sim ,xi = xi_sim, quantile =quantile, u = u_sim)
plot_quantile <- ggplot(simulation_data_df, aes(x = x_plot, y = y_plot, colour = type, fill = type)) +
  geom_line()+
  geom_area(alpha = 0.25, position = 'identity')+
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
  xlim(1.25,6)

#save plot
file_name <- paste('UQ--posterior',  quantile_value,'quantile_penultimate_parametric',
                   prior_sim, 'prior',
                   sigma_sim, 'sigma',
                   xi_sim, 'xi.png', sep = '-')

ggsave(here::here('06-measurement-scale-uncertainty', 'figures', file_name),
       plot = plot_quantile,
       units = "px",
       width = 1197,
       height = 683,
       dpi = 300)


