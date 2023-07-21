#############################################
#Likelihood function that parametrically account for uncertainty in the measurement scale
#############################################


#Likelihood function for the threshold exceedance model under penultimate
#approximation
llh_gpd_scale_penultimate <- function(sigxi, lambda, c, u, x, negative = FALSE){

  # Check inputs
  stopifnot(exprs = {
    is.numeric(u)
    length(u) == 1
    is.numeric(lambda)
    length(lambda) == 1
    is.numeric(c)
    length(c) == 1
    is.numeric(sigxi)
    length(sigxi)==2 #unique shape and scale_u parameters
    is.logical(negative)
    u <= min(v)
  })

  # Latent GPD parameters
  sig <- sigxi[1]
  sig_y <- sig * u ^(lambda-1)
  xi <- sigxi[2]
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

  return((-1)^negative * llh)
}


#Likelihood function under the assumption that the GP tail form applies for a
#range of choices of measurement scale lambda
llh_gpd_scale <- function(sigxi, lambda, u, x, negative = FALSE){

  # Check inputs
  stopifnot(exprs = {
    is.numeric(u)
    length(u) == 1
    is.numeric(lambda)
    length(lambda) == 1
    is.numeric(sigxi)
    length(sigxi)==2 #unique shape and scale_u parameters
    is.logical(negative)
    u <= min(v)
  })

  # Latent GPD parameters
  sig <- sigxi[1]
  xi <- sigxi[2]

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
  prob <- x^(lambda-1) * dgpd(x = x_mod, shape = xi, scale = sig, shift = u_mod)
  llh <- sum(log(prob))

  return((-1)^negative * llh)
}
