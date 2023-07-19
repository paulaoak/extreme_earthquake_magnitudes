###########################
# LOG-POSTERIOR
###########################


###########################
# CONSTANT THRESHOLD UNROUNDED OBSERVATIONS
###########################
log_posterior_gpd_mmax_mean <- function(params, x, prior, threshold, upper_mmax = NULL, b_value=NULL, epsilon=NULL, alpha1 = NULL, beta1 = NULL, alpha2 = NULL, beta2 = NULL) {
  log_posterior <- llh_gpd_mmax_mean_varu(mmax_mean = params, u = threshold, v = threshold, x = x)
  + prior(params, threshold = threshold, b_value = b_value, epsilon = epsilon, alpha1 = alpha1, beta1 = beta1, upper_mmax = upper_mmax, alpha2 = alpha2, beta2 = beta2)

  return(log_posterior)
}

