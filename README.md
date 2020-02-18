# Oxford mini-seminar on Bayesian statistics

Code associated with the introduction to Bayesian statistics held on the morning of Tuesday, February 18th 2020, at the Department of Zoology of Oxford University. 

Start from `model_mean.R`, which simulates and fits the easiest of data. This script has examples of moment matching, and of how to get rid of attenuation bias (due to uncertainty in the predictor) in a linear regression.

The `ipm_uncertainty.R` script simulates data on survival, growth, and reproduction to parameterize an IPM, builds the IPM fitting linear models on the simulated demographic data, and then estimates uncertainty in the asymptotic population growth rate ($\lambda) using bootstrap sampling, the confidence intervals of linear models, and 1000 posterior samples from a Bayesian model of the three linear models (those on survival, growth, and reproduction).
