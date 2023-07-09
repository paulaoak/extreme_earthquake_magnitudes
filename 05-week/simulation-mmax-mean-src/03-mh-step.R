##########################
# MH STEPS
##########################

# MH step with random walk proposal
mh_step_random_walk <- function(current_value, log_posterior, sd = 1){
  x <- current_value
  y <- rnorm(1, mean = x, sd = sd) #generate proposal

  log_ratio <- log_posterior(y)-log_posterior(x) #acceptance log_ratio

  #accept or reject the proposed values
  if (runif(1)<= exp(log_ratio)){
    x<-y
  }
  x
}

# MH step with random walk proposal with proposal only generating positive value
mh_step_random_walk_positive <- function(current_value, log_posterior, sd = 1){
  x <- current_value
  y <- truncnorm::rtruncnorm(1, mean = x, sd = sd, a = 0) #generate proposal using truncated normal distribution

  log_ratio <- (log_posterior(y) + log(truncnorm::dtruncnorm(x, a=0, mean = y, sd = sd)))-
    (log_posterior(x) + log(truncnorm::dtruncnorm(y, a=0, mean = x, sd = sd))) #acceptance log_ratio

  #accept or reject the proposed values
  if (log(runif(1))<= log_ratio){
    x<-y
  }
  x
}
