data{
  
  int n;
  vector[n] x;
  vector[n] y;
  
}

parameters {
  real beta0;
  real beta1;
  real<lower=0> y_sd;
}

transformed parameters{
  vector[n] yhat;
  
  yhat = beta0 + beta1 * x;
}

model {
  
  // priors
  beta0 ~ normal(0,1);
  beta1 ~ normal(0,1);
  y_sd  ~ inv_gamma(0.1, 0.1);
  
  // model
  y ~ normal(yhat, y_sd);
  
}
