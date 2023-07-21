#############################################
#MLE when parametrically accounting for uncertainty in the measurement scale
#############################################


#MLE for the threshold exceedance model under penultimate approximation
mle_gpd_scale_penultimate <- function(sigxi_lambda_c, u, x, llh_val = TRUE, hessian = FALSE, method =  "Nelder-Mead", maxit = 10000, ...){

  # Check inputs
  stopifnot(exprs = {
    is.numeric(u)
    length(u) == 1
    is.numeric(sigxi_lambda_c)
    length(sigxi_lambda_c)==4 #unique shape, scale, lambda and c parameters
    is.logical(llh_val)
    is.logical(hessian)
  }
  )

  # GPD parameters
  sig_0 <- sigxi[1]
  xi_0 <- sigxi[2]

  # Check valid starting point
  condition <- (xi_0 < 0) && (max(x) >= (u - sig_0 / xi_0))
  if(condition){stop('Invalid starting point. Data above upper end point of distribution.')}
  stopifnot(sig_0>0)

  # Numerically minimise negative log-likelihood
  estimate <- optim(fn = llh_gpd_scale_penultimate,
                    par = sigxi_lambda_c,
                    u = u,
                    x = x,
                    negative = TRUE,
                    hessian = hessian,
                    method = method,
                    control = list(maxit = maxit, ...))

  if(estimate$convergence)
    warning("optimization may not have succeeded")

  # Format output
  if(!llh_val & !hessian){out <-  estimate$par} else {out <- list(params = estimate$par)}
  if(llh_val) out$loglik <- -estimate$value
  if(hessian) out$hessian <- estimate$hessian

  return(out)
}


#Likelihood function under the assumption that the GP tail form applies for a
#range of choices of measurement scale lambda
mle_gpd_scale_penultimate <- function(sigxi_lambda, u, x, llh_val = TRUE, hessian = FALSE, method =  "Nelder-Mead", maxit = 10000, ...){

  # Check inputs
  stopifnot(exprs = {
    is.numeric(u)
    length(u) == 1
    is.numeric(sigxi_lambda)
    length(sigxi_lambda)==3 #unique shape, scale and lambda parameters
    is.logical(llh_val)
    is.logical(hessian)
  }
  )

  # GPD parameters
  sig_0 <- sigxi[1]
  xi_0 <- sigxi[2]

  # Check valid starting point
  condition <- (xi_0 < 0) && (max(x) >= (u - sig_0 / xi_0))
  if(condition){stop('Invalid starting point. Data above upper end point of distribution.')}
  stopifnot(sig_0>0)

  # Numerically minimise negative log-likelihood
  estimate <- optim(fn = llh_gpd_scale,
                    par = sigxi_lambda,
                    u = u,
                    x = x,
                    negative = TRUE,
                    hessian = hessian,
                    method = method,
                    control = list(maxit = maxit, ...))

  if(estimate$convergence)
    warning("optimization may not have succeeded")

  # Format output
  if(!llh_val & !hessian){out <-  estimate$par} else {out <- list(params = estimate$par)}
  if(llh_val) out$loglik <- -estimate$value
  if(hessian) out$hessian <- estimate$hessian

  return(out)
}
