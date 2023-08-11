#######################
##Tapered GPD with constant threshold Bayesian estimation for uninformative uniform prior
#######################

# LOG-PRIOR---------------------------------------------------------------------
# Define the uninformative uniform prior function for the tapered GPD parameters
unif_log_prior <- function(params) {
  beta <- params[1]
  theta <- params[2]

  if(beta<0 || theta<0){
    log_prior <- -1e7
  }
  else{
    #uniform prior for log(beta) and log(theta)
    log_prior <- - log(beta) - log(theta)
  }

  return(log_prior)
}
