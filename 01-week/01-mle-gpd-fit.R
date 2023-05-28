##############################
## Maximum likelihood fitting for GPD with constant (predefined threshold)
##############################

# Required libraries
library(ismev)

# Selected threshold
u <- 30

fit_gpd <- function(data, threshold, method =  "Nelder-Mead", maxit = 10000, ...){
  
  #compute excesses
  excess <- data[data>threshold]-threshold
  
  #number of observations exceeding the threshold
  Nu <- length(excess)
  
  #log likelihood(parameters); sigma = parameter[1] and xi = parameter[2]
  minus_likelihood <- function(parameter){
    sigma <- parameter[1]
    xi <- parameter[2]
    
    #satisfy conditions to be a valid GPD distribution
    condition_1 <- sigma <= 0
    condition_2 <- (xi < 0) && (max(excess) > ( - sigma/xi))
    
    if(condition_1 || condition_2)
      minus_log_lik <- 1e6
    else {
      minus_log_lik <- Nu*log(sigma) + (1+1/xi)*sum(log(1+xi*excess/sigma))
    }
    minus_log_lik
}
  
  #initial mle estimate (based on the true mean and variance of GPD distribution)
  x_mean <- mean(excess)
  x_var <- var(excess)
  xi_0 <- -0.5 * (((x_mean * x_mean)/x_var) - 1)
  sigma_0 <- 0.5 * x_mean * (((x_mean * x_mean)/x_var) + 1)
  initial_estimate <- c(sigma_0, xi_0) 
  
  #we use the optim function to maximise the log likelihood
  mle_estimates <- optim(initial_estimate, 
                         minus_likelihood, #recall that optim minimises
                         hessian = TRUE,
                         method = method,
                         control = list(maxit = maxit, ...))
  if(mle_estimates$convergence)
    warning("optimization may not have succeeded")
  
  #compute inverse fisher information matrix to obtain standard error of mles
  inverse_fisher <- solve(mle_estimates$hessian)
  
  return(list(mle = mle_estimates$par,
              standard_errors = sqrt(diag(inverse_fisher)),
              maximum_log_likelihood = -mle_estimates$value)) 
}

  
#example
data(rain)
fit_gpd(rain, u)
