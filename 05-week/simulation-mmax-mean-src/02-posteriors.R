###########################
# LOG-POSTERIOR
###########################


###########################
# CONSTANT THRESHOLD UNROUNDED OBSERVATIONS
###########################
log_posterior_gpd_mmax_mean <- function(params, x, prior, threshold, upper_mmax = NULL, b_value, epsilon, alpha = NULL, beta = NULL) {
  log_posterior <- llh_gpd_mmax_mean_varu(mmax_mean = params, u = threshold, v = threshold, x = x)
  + prior(params, threshold = threshold, b_value = b_value, epsilon = epsilon, alpha = alpha, beta = beta, upper_mmax = upper_mmax)

  return(log_posterior)
}

