##############################
## Fit a GPD with constant (predefined) threshold using Bayesian estimation
##############################

library(VGAM)
library(MASS)

## Example 1-------------------------------------------------------------
# Generate sample data
set.seed(123)
data <- rgpd(100, shape = 0.1, scale = 1)

# Define the log-likelihood function for the GPD
log_likelihood <- function(params, data, threshold = 0) {
  xi <- params[1]
  sig_u <- params[2]
  Nu <- length(data)
  excess<- data - threshold
  
  log_likelihood <- -Nu*log(sig_u) - (1+1/xi)*sum(log(1+xi*excess/sig_u))
  
  return(log_likelihood)
}

# Define the log-prior function for the GPD parameters
log_prior <- function(params) {
  xi <- params[1]
  sig_u <- params[2]
  
  #uniform priors for both parameters
  log_prior <- dunif(xi, min = 0, max = 2, log = TRUE) + 
    dunif(sig_u, min = 0, max = 3, log = TRUE)
  
  return(log_prior)
}

# Define the log-posterior function for the GPD
log_posterior <- function(params, data, threshold = 0) {
  log_posterior <- log_likelihood(params, data, threshold) + log_prior(params)
  
  return(log_posterior)
}

# Perform MCMC sampling using the Metropolis-Hastings algorithm
n_iter <- 1e5
n_burn <- 1000
n_chains <- 4

# Initial parameter values
params_current <- c(xi = 0.15, sig_u = 0.75)

# Storage for sampled parameter values
params_samples <- matrix(NA, nrow = (n_iter - n_burn) * n_chains, ncol = 2)

# Metropolis-Hastings algorithm
for (chain in 1:n_chains) {
  for (iter in 1:n_iter) {
    #generate proposed parameter values from a multivariate normal distribution
    params_proposed <- mvrnorm(1, params_current, matrix(c(0.11, -0.02, -0.02, 0.02), nrow = 2))
    params_proposed[1] <- max(params_proposed[1], 0.00001)
    params_proposed[2] <- max(params_proposed[2], 0.00001)
    
    #acceptance ratio
    log_ratio <- log_posterior(params_proposed, data, 0) - 
      log_posterior(params_current, data, 0)
    ratio <- exp(log_ratio)
    
    #accept or reject the proposed values
    if (runif(1) < ratio) {
      params_current <- params_proposed
    }
    
    #store the parameter values after burn-in
    if (iter > n_burn) {
      params_samples[(iter - n_burn) + (chain - 1) * (n_iter - n_burn), ] <- params_current
    }
  }
}

# Summary of the sampled parameters
params_summary <- apply(params_samples, 2, function(x) {
  c(mean = mean(x), sd = sd(x), quantile(x, c(0.025, 0.975)))
})
print(params_summary)

# Check autocorrelation functions
acf(params_samples[,1])
acf(params_samples[,2])
