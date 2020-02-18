data{
  int n;
  vector[n] x;
}

parameters{
  real u;
  real<lower=0> y_sd;
}

// transformed parameters{
//   real A;
//   real B;
// 
//   A = square(u) / y_sd;
//   B = u / y_sd;
// 
// }

model{
  
  // priors
  u ~ normal(0, 100);
  y_sd ~ gamma(0.01, 0.01);
  
  // sampling
  x ~ gamma( square(u) / y_sd, u / y_sd);
}
