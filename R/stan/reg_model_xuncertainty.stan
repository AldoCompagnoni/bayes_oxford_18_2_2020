data{
  
  int n;
  vector[n] x1;
  vector[n] x2;
  vector[n] x3;
  vector[n] x4;
  vector[n] x5;
  vector[n] x6;
  vector[n] x7;
  vector[n] x8;
  vector[n] x9;
  vector[n] x10;
  vector[n] y;
  
}

parameters {
  real y_int;
  real y_b;
  real<lower=0> y_sd;
  real<lower=0> x_sd;
  vector[n] x_true;
}

transformed parameters{
  vector[n] yhat;
  
  yhat = y_int + y_b * x_true;
}

model {
  
  // placeholder
  int id_i;
  
  // priors
  y_int ~ normal(0,1);
  y_b   ~ normal(0,1);
  y_sd  ~ inv_gamma(0.1, 0.1);
  x_sd  ~ inv_gamma(0.1, 0.1);
  
  // models
  x1 ~ normal(x_true, x_sd);
  x2 ~ normal(x_true, x_sd);
  x3 ~ normal(x_true, x_sd);
  x4 ~ normal(x_true, x_sd);
  x5 ~ normal(x_true, x_sd);
  x6 ~ normal(x_true, x_sd);
  x7 ~ normal(x_true, x_sd);
  x8 ~ normal(x_true, x_sd);
  x9 ~ normal(x_true, x_sd);
  x10 ~ normal(x_true, x_sd);
  y ~ normal(yhat, y_sd);
  
}
