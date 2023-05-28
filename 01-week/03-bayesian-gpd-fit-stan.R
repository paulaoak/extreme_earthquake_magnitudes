library(rstan)

# Sample data
library(ismev)
data(rain)
data <- rain[rain>30]

# Define the Stan model
stan_code <- "
data {
  int<lower=0> N;
  real<lower=0> threshold;
  real<lower=0> y[N];
}
parameters {
  real xi;
  real<lower=0> sig_u;
}
model {
  xi ~ uniform(-1, 1);
  sig_u ~ uniform(0.001, 13);
  
  for (n in 1:N) {
    target += - log(sig_u) - (1 + 1/xi) * log((1+ xi*(y[n] - threshold)/sig_u));
  }
}
"

# Prepare the data
stan_data <- list(N = length(data),
                  threshold = 30,
                  y = data)

# Compile the model
model <- stan_model(model_code = stan_code)

# Fit the model using MCMC sampling
fit <- sampling(model, data = stan_data, iter = 5e4, warmup = 1000, chains = 4, cores = 4) #parallelise computations using several cores

# Summary of the fitted parameters
print(fit)

# Diagnostic plots
plot(fit)

traceplot(fit, pars = c("xi", "sig_u"), inc_warmup = TRUE, nrow = 2)

# Further inspection of all chains combined
sampler_params <- get_sampler_params(fit, inc_warmup = TRUE)
summary(do.call(rbind, sampler_params), digits = 2)
