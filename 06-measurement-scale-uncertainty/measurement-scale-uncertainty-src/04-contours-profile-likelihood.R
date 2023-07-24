xi <- -0.084
sigma <- 0.48
u <- 1
n_data <- 500
x <- rgpd(n = n_data, scale = sigma, shape = xi, shift = u)
x_mod <- log(x)
mle_estimates <- mle_gpd_scale_penultimate_profile(sigxi_lambda = c(0.48,-0.084, 0.1), u = 1, x = x)

llh_gpd_scale_penultimate_for_profile(sigxi_lambda = c(0.52960778, -0.01042333, 1.50886393), u = u, x = x)

llh_gpd_scale_penultimate_for_profile(sigxi_lambda = c(0.48, -0.084, 1), u = u, x = x)

library(ggplot2)
library(dplyr)

input_1 <- seq(0, 3,0.05)
input_2 <- seq(-0.35, 0.35,0.01)

df <- expand.grid(input_1,input_2)
llh_value <- numeric(dim(df)[1])
for (i in 1:dim(df)[1]){
  lambda_val <- df$Var1[i]
  xi_val <- df$Var2[i]
  llh_value[i] <- profile_mle_xi_lambda_new(lambda = lambda_val, xi = xi_val, u = u, x = x, sig = 0.1, method = 'Brent')
}

df$z <- llh_value
ggplot(df, aes(x = Var1, y = Var2, z=z)) +
  geom_contour_filled()


library(plotly)

input_1 <- seq(-1, 3.4,0.05)
input_2 <- seq(-0.15, 1.95, 0.01)

llh_value <- matrix(data = 0, nrow = length(input_1), ncol = length(input_2))
for(i in 1:length(input_1)){
  for(j in 1:length(input_1)){
    llh_value[i,j] <- profile_mle_xi_lambda_new(lambda = input_1[i], xi = input_2[j], u = u, x = x, sig = 0.1, method = 'Brent')
  }
}

plot_ly(x = input_1, y = input_2, z = llh_value) %>% add_surface()



my_function <- function(x,y) {

  final_value = x^2 + y^2
}


input_1 <- seq(-1.5, 1.7,0.1)
input_2 <- seq(-1.5, 1.5,0.1)


z <- outer(input_1, input_2, my_function)


contours_profile_llh <- contourlevel(
  profile_mle_xi_lambda,
  p = c(0.6826895, 0.9544997),
  xmin = c(-1, -1),
  xmax = c(3.5, 0.95),
  neval = 10000,
  napprox = 10,
  u = u,
  x = x,
  sig = 0.45
)



lambda_vec <- seq(0, 3, 0.05)
xi_vec <- seq(-0.35, 0.35, 0.01)
df <- expand.grid(lambda_vec = lambda_vec, xi_vec = xi_vec)

#define weights to use
llh_value <- numeric(dim(df)[1])
for (i in 1:dim(df)[1]){
  lambda_val <- df$lambda_vec[i]
  xi_val <- df$xi_vec[i]
  llh_value[i] <- profile_mle_xi_lambda_new(lambda = lambda_val, xi = xi_val, u = u, x = x, sig = 0.1, method = 'Brent')
}

wt <- exp(-2*(mle_estimates$loglik + llh_value))

df$lambda_vec_mod <- df$lambda_vec-1
#perform weighted least squares regression
wls_model <- lm(xi_vec ~ lambda_vec_mod, data = df, weights=wt)

