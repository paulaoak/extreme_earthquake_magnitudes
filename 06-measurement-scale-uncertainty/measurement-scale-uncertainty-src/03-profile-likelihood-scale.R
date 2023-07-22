#############################################
#MLE when parametrically accounting for uncertainty in the measurement scale
#############################################


#MLE for the threshold exceedance model under penultimate approximation
profile_mle_c <- function(u, x, sigxi, lambda_min_grid, lambda_max_grid, step_grid, llh_val = TRUE, hessian = FALSE, method =  "Nelder-Mead", maxit = 10000, ...){

  # Check inputs
  stopifnot(exprs = {
    is.numeric(u)
    length(u) == 1
    is.numeric(lambda_min_grid)
    length(lambda_min_grid) == 1
    is.numeric(lambda_max_grid)
    length(lambda_max_grid) == 1
    is.numeric(step_grid)
    length(step_grid) == 1
    is.numeric(sigxi)
    length(sigxi)==2
    is.logical(llh_val)
    is.logical(hessian)
  }
  )

  # Create the grid of lambda values
  lambda_grid <- seq(lambda_min_grid, lambda_max_grid, by = step_grid)

  # GPD parameters
  sig_0 <- sigxi[1]
  xi_0 <- sigxi[2]

  # Check valid starting point
  condition <- (xi_0 < 0) && (max(x) >= (u - sig_0 / xi_0))
  if(condition){stop('Invalid starting point. Data above upper end point of distribution.')}
  stopifnot(sig_0>0)

  # Pre-allocate space to vector with profile MLE estimations of xi_y
  out <- numeric(length(lambda_grid))

  for(i in 1:length(lambda_grid)){
    # Numerically minimise negative profile log-likelihood for fixed lambda
    estimate <- optim(fn = profile_llh_gpd_scale_penultimate,
                      par = sigxi,
                      u = u,
                      x = x,
                      lambda = lambda,
                      negative = TRUE,
                      hessian = hessian,
                      method = method,
                      control = list(maxit = maxit, ...))
    if(estimate$convergence)
      warning("optimization may not have succeeded")
    # Format output
    out <-  estimate$par[2]
  }
  return(out)
}


# Estimate c using profile log-likelihood
#MLE for the threshold exceedance model under penultimate approximation
profile_mle_xi_lambda <- function(lambda_xi, u, x, sig, method =  "Nelder-Mead", maxit = 10000, ...){

  # Check inputs
  stopifnot(exprs = {
    is.numeric(u)
    length(u) == 1
    is.numeric(lambda_xi)
    length(lambda_xi)==2
  }
  )

  # GPD parameters
  xi_y <- lambda_xi[2]
  lambda <- lambda_xi[1]
  sig_y_0 <- sig * u^(lambda-1)

  if(lambda==0){
    x_mod <- log(x)
    u_mod <- log(u)
  }
  else{
    x_mod <- (x^(lambda)-1)/lambda
    u_mod <- (u^(lambda)-1)/lambda
  }

  # Check valid starting point
  condition <- (xi_y < 0) && (max(x_mod) >= (u_mod - sig_y_0 / xi_y))
  if(condition){stop('Invalid starting point. Data above upper end point of distribution.')}
  stopifnot(sig_y_0>0)

  # Numerically minimise negative profile log-likelihood for fixed lambda
  estimate <- optim(fn = profile_llh_gpd_scale_penultimate,
                    par = sig,
                    u = u,
                    x = x,
                    xi = xi_y,
                    lambda = lambda,
                    negative = TRUE,
                    hessian = FALSE,
                    method = method,
                    control = list(maxit = maxit, ...))
  if(estimate$convergence)
    warning("optimization may not have succeeded")

  # Format output
  out <-  estimate$value

  return(out)
}


profile_mle_xi_lambda_new <- function(lambda, xi, u, x, sig, method =  "Nelder-Mead", maxit = 10000, ...){

  # Check inputs
  stopifnot(exprs = {
    is.numeric(u)
    length(u) == 1
    is.numeric(xi)
    length(xi)==1
    is.numeric(lambda)
    length(lambda)==1
  }
  )

  # GPD parameters
  xi_y <- xi
  sig_y_0 <- sig * u^(lambda-1)

  if(lambda==0){
    x_mod <- log(x)
    u_mod <- log(u)
  }
  else{
    x_mod <- (x^(lambda)-1)/lambda
    u_mod <- (u^(lambda)-1)/lambda
  }

  # Check valid starting point
  condition <- (xi_y < 0) && (max(x_mod) >= (u_mod - sig_y_0 / xi_y))
  #if(condition){stop('Invalid starting point. Data above upper end point of distribution.')}
  stopifnot(sig_y_0>0)

  # Numerically minimise negative profile log-likelihood for fixed lambda
  estimate <- optim(fn = profile_llh_gpd_scale_penultimate,
                    par = sig,
                    u = u,
                    x = x,
                    xi = xi_y,
                    lambda = lambda,
                    negative = TRUE,
                    hessian = FALSE,
                    method = method,
                    lower = 0.01,
                    upper = 5,
                    control = list(maxit = maxit, ...))
  if(estimate$convergence)
    warning("optimization may not have succeeded")

  # Format output
  out <-  estimate$value

  return(out)
}

