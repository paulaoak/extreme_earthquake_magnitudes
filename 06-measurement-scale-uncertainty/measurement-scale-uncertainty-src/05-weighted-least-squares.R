##################################
## Estimate parameter c of GPD model accounting for measurement scale uncertainty
## under penultimate approximation
##################################

estimate_c_penultimate <- function(min_lambda, max_lambda, step_lambda,
                                   min_xi, max_xi, step_xi, u, x){
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





