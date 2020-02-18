data{
  
  int n;
  int y[n];
  
}

parameters {
  real beta0;
}

transformed parameters{
  real yhat;
  yhat = inv_logit(beta0);
}

model {

  // model
  y ~ bernoulli(yhat);
  
}
