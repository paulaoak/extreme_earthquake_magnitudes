#######################
##Truncated GR Bayesian estimation for uninformative uniform prior
#######################

# LOG-PRIOR---------------------------------------------------------------------
# Define the uninformative uniform prior function for the tapered GPD parameters
unif_log_prior_gr <- function(params, u) {
  beta <- params[1]
  mmax <- params[2]

  if(beta<0 || mmax<u){
    log_prior <- -1e7
  }
  else{
    #uniform prior for beta and mmax
    log_prior <- 0
  }

  return(log_prior)
}
