# Simulate simple data and fit using Bayesian models in STAN
# Sections are numbered
# 1. Normal distribution
# 2. Normal regression and attenuation bias
# 3. Gamma model: moment matching
library(rstan)
library(dplyr)
library(tidyr)
library(shinystan)

# set rstan options
rstan_options( auto_write = TRUE )
options( mc.cores = parallel::detectCores() )

# 1. Normal distribution -----------------------------------------
x <- rnorm( 1000, 5 )

data_stan <- list(
  n = length(x),
  x = x
)

# simulation parameters
sim_pars <- list(
  warmup = 1000, 
  iter = 4000, 
  thin = 2,
  chains = 3 
)

# fit the normal distribution
fit  <- stan(
  file = 'R/stan/mean_model.stan',
  data = data_stan,
  pars = c('u','y_sd'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains
)

shinystan::launch_shinystan(fit)


# 2. Normal regression and attenuation bias ------------------------

# simulate 
x  <- rnorm( 1000, 5 )
y  <- rnorm( 1000, 2 + x*1, 0.5 )

# play with the uncertainty in the detection of x
# IMPOARTANT: this parameter is SET TO 1 in the stan model
# hence, if you change this in here, you'll need to either 
bias_sd <- 0.5

# produce data, simulating impefect detection of 
data_stan <- list(
  n  = length(x),
  # this is a 
  x  = rnorm( 1000, x, bias_sd ),
  x1 = rnorm( 1000, x, bias_sd ),
  x2 = rnorm( 1000, x, bias_sd ),
  x3 = rnorm( 1000, x, bias_sd ),
  x4 = rnorm( 1000, x, bias_sd ),
  x5 = rnorm( 1000, x, bias_sd ),
  x6 = rnorm( 1000, x, bias_sd ),
  x7 = rnorm( 1000, x, bias_sd ),
  x8 = rnorm( 1000, x, bias_sd ),
  x9 = rnorm( 1000, x, bias_sd ),
  x10 = rnorm( 1000, x, bias_sd ),
  y  = y
)


# correct estimates
lm(y ~ x)

# estimates modified by attenuation bias
lm(data_stan$y ~ data_stan$x1)

plot(y ~ x)
plot(data_stan$x1,data_stan$y)


# regression model with attenuation bias
fit_bias  <- stan(
  file = 'R/stan/reg_model.stan',
  data = data_stan,
  pars = c('y_int','y_b'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains
)

shinystan::launch_shinystan(fit)

# fit (hopefully unbiased) model!
fit_unbias  <- stan(
  file = 'R/stan/reg_model_xuncertainty.stan',
  data = data_stan,
  pars = c('y_int','y_b','y_sd','x_sd'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains
)


# 3. Gamma model: moment matching ----------------------------

# simulate data
x <- rgamma(1000, 0.5, 0.2)

# store info in a list
data_stan <- list(
  n = length(x),
  x = x
)

# recover parameters
# look in stan model to examine moment matching
fit  <- stan(
  file = 'R/stan/mean_gamma_model.stan',
  data = data_stan,
  pars = c('u','y_sd'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains
)

