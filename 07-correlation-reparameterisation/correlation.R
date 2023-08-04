# Import necessary libraries
library(ggplot2)
library(latex2exp)

#Function to compute correlation under different parameterisations
correlation <- function(sigma, xi){
  #expression for correlation under shape-scale parameterisation
  rho_xi_sigma <- 1/sqrt(2*(xi +1))

  #expression for correlation under mean-mmax parameterisation
  mean <- sigma / (1-xi)
  mmax <- - sigma / xi
  i_xi_xi <- 2/((2*xi+1)*(xi+1))
  i_sigma_xi <- 1/(sigma * (2*xi+1) * (xi+1))
  i_sigma_sigma <- 1/(sigma^2 * (2*xi+1))

  i_mmax_mmax <- mean^2 * (i_xi_xi -2 * mean * i_sigma_xi + mean^2 * i_sigma_sigma)
  i_mean_mean <- mmax^2 * (i_xi_xi -2 * mmax * i_sigma_xi + mmax^2 * i_sigma_sigma)
  i_mmax_mean <- mmax * mean * (-i_xi_xi + (mean + mmax) * i_sigma_xi - mmax * mean * i_sigma_sigma)

  rho_mean_mmax <- i_mmax_mean /sqrt(i_mmax_mmax * i_mean_mean)

  return(c(rho_xi_sigma, rho_mean_mmax))
}

#Test values
sigma <- seq(0.05, 3.5, by = 0.01)
xi <- seq(-0.49, -0.01, by = 0.001)

corr_value_min = matrix(NA, nrow = length(sigma), ncol = length(xi))
corr_value_new = matrix(NA, nrow = length(sigma), ncol = length(xi))
corr_value_old = matrix(NA, nrow = length(sigma), ncol = length(xi))
for(i in 1:length(sigma)){
  for(j in 1:length(xi)){
    correlation_values = correlation(sigma[i], xi[j])
    corr_value_min[i,j] = correlation_values[1]>correlation_values[2]
    corr_value_new[i,j] = correlation_values[2]
    corr_value_old[i,j] = correlation_values[1]
  }
}

# Check whether correlation for mmax-mean is always smaller than that of xi-sigma
sum(corr_value_min) == length(sigma)*length(xi)
# Print max and min values for each correlation coefficient
max(corr_value_new)
min(corr_value_new)
max(corr_value_old)
min(corr_value_old)

# Construct dataframe for plotting
df <- data.frame(correlation = c(corr_value_new, corr_value_old),
                 type = rep(c('M_max-mean', 'sigma_xi'), each = length(corr_value_new)))

corr_distribution <- ggplot(df, aes(x = correlation, colour = type, fill = type)) +
  geom_histogram(alpha = 0.4, bins = 45)+
  theme_classic()+
  xlab('Correlation coefficient')+
  theme(#axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(family = 'sans'),
        axis.title.x = element_text(family = 'sans', size = 9),
        axis.text.y = element_text(family = 'sans'),
        title = element_text(family = 'sans', size = 8),
        legend.title = element_blank(),
        legend.text = element_text(family = 'sans', size = 9),
        panel.grid.major = element_line(),
        #panel.grid.minor = element_line(),
        axis.line.y = element_line(),
        axis.ticks.y = element_line(),
        axis.line.x = element_line(),
        axis.ticks.x = element_line(),
        legend.text.align = 0)+
  #legend.direction='vertical')+
  #legend.position = 'right')+
  expand_limits(y = 0) +
  coord_cartesian(expand = FALSE, clip = "off")+
  scale_color_discrete(labels = unname(TeX(c("$(\\bar{m}, M_{\\max})$", "$(\\sigma_u, \\xi)$"))))+
  scale_fill_discrete(labels = unname(TeX(c("$(\\bar{m}, M_{\\max})$", "$(\\sigma_u, \\xi)$"))))

corr_distribution



plot(corr_value_new, corr_value_old)



mean_mmax <- c(corr_value_new)
hist(mean_mmax)
sig_xi <- c(corr_value_old)
hist(sig_xi)
