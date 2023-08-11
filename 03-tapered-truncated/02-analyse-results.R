### BISECTION ALGORITHM
bisection<- function(CDF, alpha, prec){

  if(alpha >= 1 || alpha <= 0){
    stop("Error: prob must be between 0 and 1.")
  }

  # Initialize quantile function
  #We have to choose a value in the range of the CDF. After testing
  #the CDF available in R, we assume that they are define in all the
  #domain of the real numbers (that is, R).
  x<-0
  aux<-CDF(x)-alpha

  if (aux<0){
    while(aux<0){
      y<-x
      x<-x+1
      aux<-CDF(x)-alpha
    }
  }
  else{
    while(aux>=0){
      y<-x
      x<-x-1
      aux<-CDF(x)-alpha
    }
  }

  #Once we have chosen the initial values, we initialize bisection
  #method.
  a <- min(x,y)
  b <- max(x,y)
  ga <- CDF(a)-alpha
  gb <- CDF(b)-alpha

  while( b-a > prec){
    c<-(a+b)/2
    gc <- CDF(c)-alpha
    if(gc*ga < 0){
      b <- c
      gb <- gc
    }
    else{
      a <- c
      ga <- gc
    }
  }
  (a+b)/2
}

###################################
#READ DATA AND COMPUTE QUANTILES
###################################


library(ggplot2)
#library(urbnthemes)
library(magrittr)
library(dplyr)


#function to compute quantile from xi and sigma
quantile_posterior_calculation <- function(xi, sigma, quantile, u){
  p <- 1 - quantile

  quantiles_posterior <- u + sigma/xi * (p^(-xi) - 1)
  return(quantiles_posterior)
}

#function to compute quantile from xi, sigma and scale
quantile_posterior_calculation_tapered <- function(theta, beta, quantile, u){
  solution <- numeric(length(theta))
  for(i in 1:length(theta)){
    cdf_aux <- function(x){1 - (u/x)^beta[i] * exp(-(x-u)/theta[i])}

    solution[i] <- bisection(CDF = cdf_aux, prec = 0.001, alpha = quantile)
  }

  quantiles_posterior <- solution

  return(quantiles_posterior)
}


#function to compute quantile from xi, sigma and scale
quantile_posterior_calculation_truncated <- function(beta, mmax, quantile, u){

  quantiles_posterior <- u - 1/beta * log(1 - quantile * pexp(q=mmax-u, rate = beta))

  return(quantiles_posterior)
}

############################
#Analysis and visualization 1
############################
xi_sim <- -0.084
sigma_sim <- 0.48
u_sim <- 1.45
n_data_sim <- 75

#obtain data frame with results from different simulations
simulation_data_df <- data.frame()

#READ FILES

#######################
#File 1
#######################
file_name_sim_gpd <- paste('gpd_simulation', n_data_sim, 'n_data', sep = '_')

simulation_data_1 <- read.csv(here::here('03-tapered-truncated','outputs', file_name_sim_gpd),
                            colClasses=c("NULL", NA, NA), col.names = c('','sigma', 'xi'))

simulation_data <- data.frame(type = rep('GPD', length(simulation_data_1$sigma)))
simulation_data$quantile_0.5 <- quantile_posterior_calculation(xi = simulation_data_1$xi, sigma = simulation_data_1$sigma, quantile = 0.5, u = u_sim)
simulation_data$quantile_0.75 <- quantile_posterior_calculation(xi = simulation_data_1$xi, sigma = simulation_data_1$sigma, quantile = 0.75, u = u_sim)
simulation_data$quantile_0.95 <- quantile_posterior_calculation(xi = simulation_data_1$xi, sigma = simulation_data_1$sigma, quantile = 0.95, u = u_sim)
simulation_data_df <- rbind(simulation_data_df, simulation_data)

#######################
#File 2
#######################
file_name_sim_tapered <- paste('tapered_simulation',
                              n_data_sim, 'n_data',
                              sep = '_')

simulation_data_1 <- read.csv(here::here('03-tapered-truncated','outputs', file_name_sim_tapered),
                            colClasses=c("NULL", NA, NA), col.names = c('','beta', 'theta'))

simulation_data <- data.frame(type = rep('Tapered GPD', length(simulation_data_1$beta)))
simulation_data$quantile_0.5 <- quantile_posterior_calculation_tapered(beta = simulation_data_1$beta, theta = simulation_data_1$theta, quantile = 0.5, u = u_sim)
simulation_data$quantile_0.75 <- quantile_posterior_calculation_tapered(beta = simulation_data_1$beta, theta = simulation_data_1$theta, quantile = 0.75, u = u_sim)
simulation_data$quantile_0.95 <- quantile_posterior_calculation_tapered(beta = simulation_data_1$beta, theta = simulation_data_1$theta, quantile = 0.95, u = u_sim)
simulation_data_df <- rbind(simulation_data_df, simulation_data)

#######################
#File 3
#######################
file_name_sim_truncated <- paste('truncated_gr_simulation',
                                         n_data_sim, 'n_data',
                                         sep = '_')

simulation_data_1 <- read.csv(here::here('03-tapered-truncated','outputs', file_name_sim_truncated),
                            colClasses=c("NULL", NA, NA), col.names = c('','beta', 'mmax'))

simulation_data <- data.frame(type = rep('Truncated GR', length(simulation_data_1$beta)))
simulation_data$quantile_0.5 <- quantile_posterior_calculation_truncated(beta = simulation_data_1$beta, mmax = simulation_data_1$mmax, quantile = 0.5, u = u_sim)
simulation_data$quantile_0.75 <- quantile_posterior_calculation_truncated(beta = simulation_data_1$beta, mmax = simulation_data_1$mmax, quantile = 0.75, u = u_sim)
simulation_data$quantile_0.95 <- quantile_posterior_calculation_truncated(beta = simulation_data_1$beta, mmax = simulation_data_1$mmax, quantile = 0.95, u = u_sim)
simulation_data_df <- rbind(simulation_data_df, simulation_data)

#Plot 0.5 quantile
quantile<- 0.5
true_quantile <- quantile_posterior_calculation(sigma = sigma_sim ,xi = xi_sim, quantile =quantile, u = u_sim)
plot_quantile <- ggplot(simulation_data_df, aes(x = quantile_0.5, colour = type, fill = type)) +
  geom_density(alpha = 0.25)+
  geom_segment(aes(x = true_quantile, y = 0, xend = true_quantile, yend = Inf, linetype = "True quantile"), color = 'black')+
  #ggtitle(paste('Posterior distribution of quantile', quantile, sep = ' '))+
  xlab(paste(quantile, 'quantile zoom', sep = ' '))+
  theme_classic()+
  theme(#axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(family = 'sans', size = 8),
    axis.text.y = element_text(family = 'sans', size = 8),
    title = element_text(family = 'sans', size = 8),
    legend.title = element_blank(),
    legend.text = element_text(family = 'sans', size = 9),
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
  xlim(1.5,2.75)


#save plot
file_name <- paste('posterior_0.5_quantile.png', sep = '-')

ggsave(here::here('03-tapered-truncated', 'figures', file_name),
       plot = plot_quantile,
       units = "px",
       width = 1197,
       height = 683,
       dpi = 300)

#Plot 0.75 quantile
quantile<- 0.75
true_quantile <- quantile_posterior_calculation(sigma = sigma_sim ,xi = xi_sim, quantile =quantile, u = u_sim)
plot_quantile <- ggplot(simulation_data_df, aes(x = quantile_0.75, colour = type, fill = type)) +
  geom_density(alpha = 0.25)+
  geom_segment(aes(x = true_quantile, y = 0, xend = true_quantile, yend = Inf, linetype = "True quantile"), color = 'black')+
  #ggtitle(paste('Posterior distribution of quantile', quantile, sep = ' '))+
  xlab(paste(quantile, 'quantile', sep = ' '))+
  theme_classic()+
  theme(#axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(family = 'sans', size = 8),
    axis.text.y = element_text(family = 'sans', size = 8),
    title = element_text(family = 'sans', size = 8),
    legend.title = element_blank(),
    legend.text = element_text(family = 'sans', size = 9),
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
  xlim(1.7,2.75)



#save plot
file_name <- paste('posterior_0.75_quantile.png', sep = '-')

ggsave(here::here('03-tapered-truncated', 'figures', file_name),
       plot = plot_quantile,
       units = "px",
       width = 1197,
       height = 683,
       dpi = 300)


#Plot 0.95 quantile
quantile<- 0.95
true_quantile <- quantile_posterior_calculation(sigma = sigma_sim ,xi = xi_sim, quantile =quantile, u = u_sim)
plot_quantile <- ggplot(simulation_data_df, aes(x = quantile_0.95, colour = type, fill = type)) +
  geom_density(alpha = 0.25)+
  geom_segment(aes(x = true_quantile, y = 0, xend = true_quantile, yend = Inf, linetype = "True quantile"), color = 'black')+
  #ggtitle(paste('Posterior distribution of quantile', quantile, sep = ' '))+
  xlab(paste(quantile, 'quantile', sep = ' ')) +
  theme_classic()+
  theme(#axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(family = 'sans', size = 8),
    axis.text.y = element_text(family = 'sans', size = 8),
    title = element_text(family = 'sans', size = 8),
    legend.title = element_blank(),
    legend.text = element_text(family = 'sans', size = 9),
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
  xlim(2, 4.5)


#save plot
file_name <- paste('posterior_0.95_quantile.png', sep = '-')

ggsave(here::here('03-tapered-truncated', 'figures', file_name),
       plot = plot_quantile,
       units = "px",
       width = 1197,
       height = 683,
       dpi = 300)


true_0.5 <- quantile_posterior_calculation(sigma = sigma_sim ,xi = xi_sim, quantile =0.5, u = u_sim)
true_0.75 <- quantile_posterior_calculation(sigma = sigma_sim ,xi = xi_sim, quantile =0.75, u = u_sim)
true_0.95 <- quantile_posterior_calculation(sigma = sigma_sim ,xi = xi_sim, quantile =0.95, u = u_sim)

summary_error <- simulation_data_df%>%
  group_by(type)%>%
  summarise(avg_0.5 = mean(quantile_0.5),
            bias2_0.5 = (mean(quantile_0.5)-true_0.5)^2,
            var_0.5 = var(quantile_0.5),
            avg_0.75 = mean(quantile_0.75),
            bias2_0.75 = (mean(quantile_0.75)-true_0.75)^2,
            var_0.75 = var(quantile_0.75),
            avg_0.95 = mean(quantile_0.95),
            bias2_0.95 = (mean(quantile_0.95)-true_0.95)^2,
            var_0.95 = var(quantile_0.95))

summary_error <- summary_error %>%
  mutate(mse_0.5 = var_0.5 + bias2_0.5,
         mse_0.75 = var_0.75 + bias2_0.75,
         mse_0.95 = var_0.95 + bias2_0.95)
