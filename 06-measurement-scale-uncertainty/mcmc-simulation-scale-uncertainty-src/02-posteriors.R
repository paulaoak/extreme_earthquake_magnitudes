###########################
# LOG-POSTERIOR FOR GPD MODEL ACCOUNTING FOR MEASUREMENT SCALE UNCERTAINTY
###########################

###########################
# UNDER PENULTIMATE APPROXIMATION
###########################
log_posterior_gpd <- function(params, c, x, u, prior) {
  new_params <- c(params, c)
  log_posterior <- llh_gpd_scale_penultimate(sigxi_lambda_c = new_params, u = u, x = x)
  + prior(params)

  return(log_posterior)
}


###########################
# UNDER THE ASSUMPTION THAT GPD MODEL IS VALID FOR DIFFERENT LAMBDAS
###########################
log_posterior_gpd <- function(params, x, u, prior) {
  log_posterior <- llh_gpd_scale(sigxi_lambda = params, u = u, x = x)
  + prior(params)

  return(log_posterior)
}
