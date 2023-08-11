###########################
# CONSTANT THRESHOLD UNROUNDED OBSERVATIONS NO SCALE
###########################
log_posterior_gpd <- function(params, x, u, prior) {
  log_posterior <- llh_gpd_varu(sigxi = params, u = u, v = u, x = x)
  + prior(params)

  return(log_posterior)
}
