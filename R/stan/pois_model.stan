data{
  
  int n;
  vector[n] x;
  int y[n];
  
}

parameters {
  real beta0;
  real beta1;
}

transformed parameters{
  vector[n] yhat;
  
  yhat = exp(beta0 + beta1 * x);
}

model {
  
  // priors
  beta0 ~ normal(0,1);
  beta1 ~ normal(0,1);
  
  // model
  y ~ poisson(yhat);
  
}
