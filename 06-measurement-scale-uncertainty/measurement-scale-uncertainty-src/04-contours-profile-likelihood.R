

contours_profile_llh <- contourlevel(
  profile_mle_xi_lambda,
  p = c(0.6826895, 0.9544997),
  xmin = NULL,
  xmax = NULL,
  neval = 10000,
  napprox = 10,
  u = u,
  x = x,
  sig = sig_0,
)
