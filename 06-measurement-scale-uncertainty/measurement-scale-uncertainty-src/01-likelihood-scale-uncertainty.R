#############################################
#Likelihood function that parametrically account for uncertainty in the measurement scale
#############################################


#Likelihood function for the threshold exceedance model under penultimate
#approximation
llh_gpd_scale_penultimate <- function(sigxi_lambda_c, u, x, negative = FALSE){

  # Check inputs
  stopifnot(exprs = {
    is.numeric(u)
    length(u) == 1
    is.numeric(sigxi_lambda_c)
    length(sigxi_lambda_c)==4 #unique shape, scale_u, lambda and c parameters
    is.logical(negative)
  })

  # GPD parameters
  lambda <- sigxi_lambda_c[3]
  c <- sigxi_lambda_c[4]
  sig <- sigxi_lambda_c[1]
  sig_y <- sig * u ^(lambda-1)
  xi <- sigxi_lambda_c[2]
  xi_y <- xi + c*(lambda-1)


  if(lambda==0){
    x_mod <- log(x)
    u_mod <- log(u)
  }
  else{
    x_mod <- (x^(lambda)-1)/lambda
    u_mod <- (u^(lambda)-1)/lambda
  }

  # Check x all be above u
  lep_fail <- any(x < u)
  if(lep_fail){stop('Lower endpoint failure:  !all(x > u).')}

  # Function body

  # Satisfy conditions to be a valid GPD distribution:
  #check all scale parameters are positive
  condition_1 <- sig <= 0
  #check all x below upper end-point (UEP) (if it exists)
  condition_2 <- (xi < 0) && (max(x) >= (u - sig / xi))

  if(condition_1 || condition_2){
    llh <- -10e6
    return((-1)^negative * llh)}

  # Log-likelihood
  prob <- x^(lambda-1) * dgpd(x = x_mod, shape = xi_y, scale = sig_y, shift = u_mod)
  llh <- sum(log(prob))
  if(prod(prob) == 0){
    llh <- -10e6
  }

  return((-1)^negative * llh)
}


#Profile likelihood function given lambda for the threshold exceedance model under
#penultimate approximation
profile_llh_gpd_scale_penultimate <- function(sig, xi, lambda, u, x, negative = FALSE){

  # Check inputs
  stopifnot(exprs = {
    is.numeric(u)
    length(u) == 1
    is.numeric(lambda)
    length(lambda) == 1
    is.numeric(sig)
    length(sig)==1
    is.numeric(xi)
    length(xi)==1
    is.logical(negative)
  })

  # GPD parameters
  sig_y <- sig * u ^(lambda-1)
  xi_y <- xi

  if(lambda==0){
    x_mod <- log(x)
    u_mod <- log(u)
  }
  else{
    x_mod <- (x^(lambda)-1)/lambda
    u_mod <- (u^(lambda)-1)/lambda
  }

  # Check x all be above u
  lep_fail <- any(x < u)
  if(lep_fail){stop('Lower endpoint failure:  !all(x > u).')}

  # Function body

  # Satisfy conditions to be a valid GPD distribution:
  #check all scale parameters are positive
  condition_1 <- sig <= 0
  #check all x below upper end-point (UEP) (if it exists)
  condition_2 <- (xi_y < 0) && (max(x_mod) >= (u_mod - sig_y / xi_y))

  if(condition_1 || condition_2){
    llh <- -10e6
    return((-1)^negative * llh)}

  # Log-likelihood
  prob <- x^(lambda-1) * dgpd(x = x_mod, shape = xi_y, scale = sig_y, shift = u_mod)
  llh <- sum(log(prob))

  return((-1)^negative * llh)
}

#Maximum likelihood function given for profile likelihood
llh_gpd_scale_penultimate_for_profile <- function(sigxi_lambda, u, x, negative = FALSE){

  # Check inputs
  stopifnot(exprs = {
    is.numeric(u)
    length(u) == 1
    is.numeric(sigxi_lambda)
    length(sigxi_lambda)==3
    is.logical(negative)
  })

  # GPD parameters
  lambda <- sigxi_lambda[3]
  sig <- sigxi_lambda[1]
  sig_y <- sig * u ^(lambda-1)
  xi_y <- sigxi_lambda[2]

  if(lambda==0){
    x_mod <- log(x)
    u_mod <- log(u)
  }
  else{
    x_mod <- (x^(lambda)-1)/lambda
    u_mod <- (u^(lambda)-1)/lambda
  }

  # Check x all be above u
  lep_fail <- any(x < u)
  if(lep_fail){stop('Lower endpoint failure:  !all(x > u).')}

  # Function body

  # Satisfy conditions to be a valid GPD distribution:
  #check all scale parameters are positive
  condition_1 <- sig <= 0
  #check all x below upper end-point (UEP) (if it exists)
  condition_2 <- (xi_y < 0) && (max(x_mod) >= (u_mod - sig_y / xi_y))

  if(condition_1 || condition_2){
    llh <- -10e6
    return((-1)^negative * llh)}

  # Log-likelihood
  prob <- x^(lambda-1) * dgpd(x = x_mod, shape = xi_y, scale = sig_y, shift = u_mod)
  llh <- sum(log(prob))

  return((-1)^negative * llh)
}

#Likelihood function under the assumption that the GP tail form applies for a
#range of choices of measurement scale lambda
llh_gpd_scale <- function(sigxi_lambda, u, x, negative = FALSE){

  # Check inputs
  stopifnot(exprs = {
    is.numeric(u)
    length(u) == 1
    is.numeric(sigxi_lambda)
    length(sigxi_lambda)==3 #unique shape and scale_u parameters
    is.logical(negative)
  })

  # GPD parameters
  lambda <- sigxi_lambda[3]
  sig <- sigxi_lambda[1]
  xi <- sigxi_lambda[2]

  if(lambda==0){
    x_mod <- log(x)
    u_mod <- log(u)
  }
  else{
    x_mod <- (x^(lambda)-1)/lambda
    u_mod <- (u^(lambda)-1)/lambda
  }

  # Check x all be above u
  lep_fail <- any(x < u)
  if(lep_fail){stop('Lower endpoint failure:  !all(x > u).')}

  # Function body

  # Satisfy conditions to be a valid GPD distribution:
  #check all scale parameters are positive
  condition_1 <- sig <= 0
  #check all x below upper end-point (UEP) (if it exists)
  condition_2 <- (xi < 0) && (max(x_mod) >= (u_mod - sig / xi))

  if(condition_1 || condition_2){
    llh <- -10e6
    return((-1)^negative * llh)}

  # Log-likelihood
  prob <- x^(lambda-1) * dgpd(x = x_mod, shape = xi, scale = sig, shift = u_mod)
  llh <- sum(log(prob))

  return((-1)^negative * llh)
}


#Likelihood function for the threshold exceedance model under penultimate
#approximation using quadratic approximation for the shape parameter
llh_gpd_scale_penultimate_quadratic <- function(sigxi_lambda_c1_c2, u, x, negative = FALSE){

  # Check inputs
  stopifnot(exprs = {
    is.numeric(u)
    length(u) == 1
    is.numeric(sigxi_lambda_c1_c2)
    length(sigxi_lambda_c1_c2)==5 #unique shape, scale_u, lambda and c parameters
    is.logical(negative)
  })

  # GPD parameters
  lambda <- sigxi_lambda_c1_c2[3]
  c1 <- sigxi_lambda_c1_c2[4]
  c2 <- sigxi_lambda_c1_c2[5]
  sig <- sigxi_lambda_c1_c2[1]
  sig_y <- sig * u ^(lambda-1)
  xi <- sigxi_lambda_c1_c2[2]
  xi_y <- xi + c1*(lambda-1)+ c2*(lambda-1)^2


  if(lambda==0){
    x_mod <- log(x)
    u_mod <- log(u)
  }
  else{
    x_mod <- (x^(lambda)-1)/lambda
    u_mod <- (u^(lambda)-1)/lambda
  }

  # Check x all be above u
  lep_fail <- any(x < u)
  if(lep_fail){stop('Lower endpoint failure:  !all(x > u).')}

  # Function body

  # Satisfy conditions to be a valid GPD distribution:
  #check all scale parameters are positive
  condition_1 <- sig <= 0
  #check all x below upper end-point (UEP) (if it exists)
  condition_2 <- (xi < 0) && (max(x) >= (u - sig / xi))

  if(condition_1 || condition_2){
    llh <- -10e6
    return((-1)^negative * llh)}

  # Log-likelihood
  prob <- x^(lambda-1) * dgpd(x = x_mod, shape = xi_y, scale = sig_y, shift = u_mod)
  llh <- sum(log(prob))
  if(prod(prob) == 0){
    llh <- -10e6
  }

  return((-1)^negative * llh)
}
