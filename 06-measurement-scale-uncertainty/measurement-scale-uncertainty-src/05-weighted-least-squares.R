##################################
## Estimate parameter c of GPD model accounting for measurement scale uncertainty
## under penultimate approximation
##################################

estimate_c_penultimate <- function(min_lambda, max_lambda, step_lambda,
                                   min_xi, max_xi, step_xi, u, x, sigxi_lambda_init){
  lambda_vec <- seq(min_lambda, max_lambda, step_lambda)
  xi_vec <- seq(min_xi, max_xi, step_xi)
  df <- expand.grid(lambda_vec = lambda_vec, xi_vec = xi_vec)

  #define weights for the weighted least square algorithm
  llh_value <- numeric(dim(df)[1])
  for (i in 1:dim(df)[1]){
    lambda_val <- df$lambda_vec[i]
    xi_val <- df$xi_vec[i]
    llh_value[i] <- profile_mle_xi_lambda_new(lambda = lambda_val, xi = xi_val, u = u, x = x, sig = 0.1, method = 'Brent')
  }

  #obtain MLE estimates
  mle_estimates <- mle_gpd_scale_penultimate_profile(sigxi_lambda = sigxi_lambda_init, u = u, x = x)
  wt <- exp(-2*(mle_estimates$loglik + llh_value)) #weights

  #obtain covariate for the regression
  df$lambda_vec_mod <- df$lambda_vec - 1

  #perform weighted least squares regression
  wls_model <- lm(xi_vec ~ lambda_vec_mod, data = df, weights=wt)

  #return covariate coefficient and standard error
  out <- numeric(2)
  out[1] <- wls_model$coefficients[[2]]
  out[2] <- sqrt(diag(vcov(wls_model)))[[2]]

  return(out)
}



##################################
## Estimate parameter c of GPD model accounting for measurement scale uncertainty
## under penultimate approximation using quadratic approximation
##################################

estimate_c_penultimate_quadratic <- function(min_lambda, max_lambda, step_lambda,
                                       min_xi, max_xi, step_xi, u, x, sigxi_lambda_init){
  lambda_vec <- seq(min_lambda, max_lambda, step_lambda)
  xi_vec <- seq(min_xi, max_xi, step_xi)
  df <- expand.grid(lambda_vec = lambda_vec, xi_vec = xi_vec)

  #define weights for the weighted least square algorithm
  llh_value <- numeric(dim(df)[1])
  for (i in 1:dim(df)[1]){
    lambda_val <- df$lambda_vec[i]
    xi_val <- df$xi_vec[i]
    llh_value[i] <- profile_mle_xi_lambda_new(lambda = lambda_val, xi = xi_val, u = u, x = x, sig = 0.1, method = 'Brent')
  }

  #obtain MLE estimates
  mle_estimates <- mle_gpd_scale_penultimate_profile(sigxi_lambda = sigxi_lambda_init, u = u, x = x)
  wt <- exp(-2*(mle_estimates$loglik + llh_value)) #weights

  #obtain covariates for the regression
  df$lambda_vec_mod <- df$lambda_vec - 1

  df$lambda_vec_quadratic <- (df$lambda_vec_mod)^2

  #perform weighted least squares regression
  wls_model <- lm(xi_vec ~ lambda_vec_mod + lambda_vec_quadratic, data = df, weights=wt)

  #return covariate coefficient and standard error
  out <- numeric(4)
  out[1] <- wls_model$coefficients[[2]]
  out[2] <- sqrt(diag(vcov(wls_model)))[[2]]
  out[3] <- wls_model$coefficients[[3]]
  out[4] <- sqrt(diag(vcov(wls_model)))[[3]]

  return(out)
}



###############first version
estimate_c_penultimate_old <- function(min_lambda, max_lambda, step_lambda,
                                       min_xi, max_xi, step_xi, u, x, sigxi_lambda_init){
  lambda_vec <- seq(min_lambda, max_lambda, step_lambda)
  xi_vec <- seq(min_xi, max_xi, step_xi)
  df <- expand.grid(lambda_vec = lambda_vec, xi_vec = xi_vec)

  #define weights for the weighted least square algorithm
  llh_value <- numeric(dim(df)[1])
  for (i in 1:dim(df)[1]){
    lambda_val <- df$lambda_vec[i]
    xi_val <- df$xi_vec[i]
    llh_value[i] <- profile_mle_xi_lambda_new(lambda = lambda_val, xi = xi_val, u = u, x = x, sig = 0.1, method = 'Brent')
  }

  #obtain MLE estimates
  mle_estimates <- mle_gpd_scale_penultimate_profile(sigxi_lambda = sigxi_lambda_init, u = u, x = x)
  wt <- exp(-2*(mle_estimates$loglik + llh_value)) #weights

  #obtain covariate for the regression
  df$lambda_vec_mod <- df$lambda_vec - 1

  #perform weighted least squares regression
  wls_model <- lm(xi_vec ~ lambda_vec_mod, data = df, weights=wt)

  #return covariate coefficient and standard error
  out <- numeric(2)
  out[1] <- wls_model$coefficients[[2]]
  out[2] <- sqrt(diag(vcov(wls_model)))[[2]]

  return(out)
}




