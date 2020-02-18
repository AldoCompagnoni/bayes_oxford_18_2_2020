data{
  
  int n;
  vector[n] x;
  vector[n] y;
  
}

parameters {
  real y_int;
  real y_b;
  real<lower=0> y_sd;
}

transformed parameters{
  vector[n] yhat;
  
  yhat = y_int + y_b * x;
}

model {
  
  // priors
  y_int ~ normal(0,1);
  y_b   ~ normal(0,1);
  y_sd  ~ inv_gamma(0.1, 0.1);
  
  // model
  y ~ normal(yhat, y_sd);
  
}
