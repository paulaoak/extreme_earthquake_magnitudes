# LOG-POSTERIOR-----------------------------------------------------------------
# Define the log-posterior function for the tapered GPD
log_posterior_tapered <- function(params, x, u, prior) {
  params_tgpd <- c(u/params[1], 1/params[1], params[2])
  log_posterior <- llh_tgpd_rd_varu(sigxitheta = params_tgpd, u = u, v = u, x = x)
  + prior(params = params)

  return(log_posterior)
}


#x <- rgpd(50)
#log_posterior_tapered(params = c(1, 1), x = x, u = 0, prior = unif_log_prior)
