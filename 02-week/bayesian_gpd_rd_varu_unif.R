#######################
##Bayesian estimation of rounded GPD with variable threshold for uninformative uniform prior
#######################
 
# LOG-PRIOR---------------------------------------------------------------------
# Define the uninformative uniform prior function for the GPD parameters
unif_log_prior <- function(params) {
  sigma <- params[1]
  
  #uniform prior for xi and log(sigma)
  log_prior <- -log(sigma)
  
  return(log_prior)
}

# LOG-POSTERIOR-----------------------------------------------------------------
# Define the log-posterior function for the GPD
log_posterior_gpd_rd_varu <- function(params, x, v, u, prior) {
  log_posterior <- llh_gpd_rd_varu(sigxi = params, u = u, v = v, x = x)
  + prior(params)
  
  return(log_posterior)
}

# MH STEPS----------------------------------------------------------------------
# MH step with random walk proposal
mh_step_random_walk <- function(current_value, log_posterior, sd = 1){
  x <- current_value
  y <- rnorm(1, mean = x, sd = sd) #generate proposal

  log_ratio <- log_posterior(y)-log_posterior(x) #acceptance log_ratio
  
  #accept or reject the proposed values
  #if (log(runif(1))<= log_ratio){
  if (runif(1)<= exp(log_ratio)){
    x<-y
  }
  x
}

# MH step with random walk proposal with proposal only generating positive value
mh_step_random_walk_positive <- function(current_value, log_posterior, sd = 1){
  x <- current_value
  y <- rnorm(1, mean = x, sd = sd) #generate proposal
  y <- max(y, 0.0001)

  log_ratio <- log_posterior(y)-log_posterior(x) #acceptance log_ratio
  
  #accept or reject the proposed values
  if (log(runif(1))<= log_ratio){
    x<-y
  }
  x
}

# MH WITHIN GIBBS ALGORITHM-----------------------------------------------------
# Perform MCMC sampling using the Metropolis-Hastings within Gibbs algorithm
bayesian_estimation_gpd_rd_varu <- function(sigxi_initial, x, v, u, n_iter, n_burn, sd_xi, sd_sig_u){
  # Initial value of parameters
  sig_u <- sigxi_initial[1]
  xi <- sigxi_initial[2]
  
  # Storage for sampled parameter values
  params_samples <- matrix(NA, ncol = n_iter, nrow = 2)
  
  # Metropolis-Hastings within Gibbs algorithm
  for (iter in 1:n_iter) {
    
    #MH for shape parameter  
    log_posterior_sig_u <- function (scale){log_posterior_gpd_rd_varu(params = c(scale, xi), x = x, v =v, u = u, prior = unif_log_prior)}
    sig_u <- mh_step_random_walk_positive (sig_u, log_posterior_sig_u, sd = sd_sig_u)
    
    #MH for shape parameter
    log_posterior_xi <- function (shape){llh_gpd_rd_varu(sigxi = c(sig_u, shape), u = u, v = v, x = x)}
    xi <- mh_step_random_walk (xi, log_posterior_xi, sd = sd_xi)
    
    params_samples[,iter] <- c(sig_u, xi)
  }
  
  params_samples[,n_burn:n_iter]
}

## SIMULATION------------------------------------------------------------------

# Define parameter values for Bayesian estimation
n_iter <- 1e5
n_burn <- 1e4
params_initial <- c(1, 0.01) #sig_u = 1, xi = 0.01

# Generate sample data
set.seed(123)
u <- 1.1
v <- rep(c(1.6, 1.2), each = 250)
x <- rgpd_rd(n = 500, scale = 2.1, shape = -0.1, shift = v, to_nearest = 0.1)

# Run MH within Gibbs
params_samples <- bayesian_estimation_gpd_rd_varu(sigxi_initial = params_initial, x = x, v = v, u = u,
                                n_iter, n_burn, sd_xi=0.1, sd_sig_u=0.2)

#Diagnostic plots and summaries
library(coda)
# Summary of the sampled parameters
params_summary <- apply(params_samples, 1, function(x) {
  c(mean = mean(x), sd = sd(x), quantile(x, c(0.025, 0.975)))
})
print(params_summary)

# Check autocorrelation functions
acf(params_samples[1,])
acf(params_samples[2,])

# Traceplots to examine mixing of the chains
plot(params_samples[1,], type="l")
plot(params_samples[2,], type="l")

# Obtain effective size
effectiveSize(t(params_samples))
