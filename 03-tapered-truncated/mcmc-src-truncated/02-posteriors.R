# LOG-POSTERIOR-----------------------------------------------------------------

#Define likelihood
llh_truncated_gr <- function(beta_mmax = params, u = u, x = x){
  # Check inputs
  stopifnot(exprs = {
    is.numeric(u)
    length(u) == 1
    is.numeric(beta_mmax)
    length(beta_mmax)==2 #unique shape and scale_u parameters
  })

  # Truncated GR parameters
  beta <- beta_mmax[1]
  mmax <- beta_mmax[2]

  # Check x all be above v
  lep_fail <- any(x < u)
  if(lep_fail){stop('Lower endpoint failure:  !all(x > u).')}

  # Function body

  # Satisfy conditions to be a valid GPD distribution:
  #check rate parameter is positive
  condition_1 <- beta <= 0
  #check all x below upper end-point (UEP) (if it exists)
  condition_2 <- max(x) >= mmax

  if(condition_1 || condition_2){
    llh <- -10e6
    return(llh)}

  # Log-likelihood
  prob <- dexp(x = x-u, rate = beta)/pexp(q = mmax-u, rate = beta)
  llh <- sum(log(prob))

  return(llh)

}


# Define the log-posterior function for the truncated GR law
log_posterior_truncated_gr <- function(params, x, u, prior) {
  log_posterior <- llh_truncated_gr(beta_mmax = params, u = u, x = x)
  + prior(params = params, u = u)

  return(log_posterior)
}
