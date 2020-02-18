# Simulate data to parameterize an IPM
# Sections are numbered
# 1. Simulate data
# 2. Fit models
# 3. IPM functions
# 4. Produce IPM
# 5. Uncertainty via bootstrap
# 6. Uncertainty via CI
# 7. Uncertainty via Bayes
library(rstan)
library(dplyr)
library(shinystan)
library(MASS)
library(ggplot2)


# 1. simulate data -----------------------------------

# set up individual random effects
Sig  <- matrix( c(1,-.5,-.5,1), 2, 2)
b0_v <- mvrnorm( 100, c(0,-2), Sig )

# generate data. Sizes are STANDARDIZED sizes
x    <- rnorm( 100 )
x1   <- rnorm( 100, b0_v[,1] + x*0.9 )
r0   <- rpois( 100, exp(b0_v[,2] + x*1 ) )
s1   <- rbinom( 100, 1, 0.5)
rcr  <- (x-2) # recruits are 2 standard deviations below parents

# demographic data
demo_df <- data.frame( x  = x,
                       x1 = x1,
                       s1 = s1,
                       r0 = r0 )

# growth data
ggplot( demo_df ) +
  geom_point( aes(x, x1) )

# recruitment data
ggplot( demo_df ) +
  geom_point( aes(x, r0) )


# 2. fit models --------------------------------------
gr_mod  <- lm( x1 ~ x, data = demo_df )
fec_mod <- glm( r0 ~ x, data = demo_df, family='poisson' )
sr_mod  <- glm( s1 ~ 1, data = demo_df, family='binomial' )
rec_mod <- lm( rcr ~ 1) # model of the mean for recruits

# IPM parameters 
pars  <- list( surv_b0 = coef(sr_mod)[1],
               
               grow_b0 = coef(gr_mod)[1],
               grow_b1 = coef(gr_mod)[2],
               grow_sd = summary(gr_mod)$sigma,
               
               fecu_b0 = coef(fec_mod)[1],
               fecu_b1 = coef(fec_mod)[2],
               
               recr_sz = coef(rec_mod),
               recr_sd = summary(rec_mod)$sigma,
               L       = min(rcr,na.rm=T),
               U       = max(c(demo_df$x,
                               demo_df$x1),
                             na.rm=T),
               mat_siz = 50
)


# 3. IPM functions -----------------------------------

gyx <- function(y,x,pars){
  # returns a *probability density distribution* for each x value
  return( dnorm(y, 
                mean = pars$grow_b0 + pars$grow_b1 * x, 
                sd   = pars$grow_sd) )
}

# Survival at size x
sx <- function() return( 0.5 )

# transition: Survival * growth
pyx<-function(y,x,pars){
  return( sx() * gyx(y,x,pars) )
}

# function to produce recruits sizes
recruits<-function(y,pars){
  dnorm( x    = y,
         mean = pars$recr_sz,
         sd   = pars$recr_sd )
}

# Production of 1-YO seeds in seed bank from x-sized moms
fx<-function(x,pars){
  exp(pars$fecu_b0 + pars$fecu_b1*x)
}

fyx<-function(y,x,pars){
  recruits(y,pars) * fx(x,pars)
}


# 4. Produce IPM -------------------------------------

# Kernel
kernel <- function(pars){
  
  n   <- pars$mat_siz
  L   <- pars$L 
  U   <- pars$U
  #these are the upper and lower integration limits
  h   <- (U-L)/n                   #Bin size
  b   <- L+c(0:n)*h                #Lower boundaries of bins 
  y   <- 0.5*(b[1:n]+b[2:(n+1)])   #Bins' midpoints
  #these are the boundary points (b) and mesh points (y)
  
  # Fertility matrix
  Fmat            <- matrix(0,n,n)
  # Transition matrix (Growth/survival)
  Tmat            <- matrix(0,n,n)
  # Growth matrix only
  Gmat            <- matrix(0,n,n)
  
  # Banked seeds go in top row
  Fmat            <- outer(y,y,fyx,pars)
  
  # Growth/survival transitions among cts sizes
  Tmat            <- outer(y,y,pyx,pars) * h
  
  # Growth matrix only (ONLY FOR CHECKS)
  Gmat            <- outer(y,y,gyx,pars) * h
  
  # Full Kernel is simply a summation of fertility and transition matrix
  k_yx          <- Fmat + Tmat     
  
  return(list(k_yx    = k_yx,
              Fmat    = Fmat,
              Tmat    = Tmat,
              Gmat    = Gmat,
              meshpts = y))
  
}

Re(eigen(kernel(pars)$k_yx)$value[1])


# 5. Uncertainty via bootstrap ------------------------------------------

boot_lam <- function(ii){
  
  ids       <- base::sample(1:nrow(demo_df), 
                            replace = T )
  
  demo_boot <- demo_df[ids,]
  
  gr_mod  <- lm( x1 ~ x, data = demo_boot )
  fec_mod <- glm( r0 ~ x, data = demo_boot, family='poisson' )
  sr_mod  <- glm( s1 ~ 1, data = demo_boot, family='binomial' )
  rec_mod <- lm( rcr ~ 1) # model of the mean for recruits
  
  # IPM parameters 
  pars  <- list( surv_b0 = coef(sr_mod)[1],
                 
                 grow_b0 = coef(gr_mod)[1],
                 grow_b1 = coef(gr_mod)[2],
                 grow_sd = summary(gr_mod)$sigma,
                 
                 fecu_b0 = coef(fec_mod)[1],
                 fecu_b1 = coef(fec_mod)[2],
                 
                 recr_sz = coef(rec_mod),
                 recr_sd = summary(rec_mod)$sigma,
                 L       = min(rcr,na.rm=T),
                 U       = max(c(demo_df$x,
                                 demo_df$x1),
                               na.rm=T),
                 mat_siz = 50
  )
  
  Re(eigen(kernel(pars)$k_yx)$value[1])
  
}
  
lam_boot_v <- sapply(1:1000, boot_lam )

boxplot(lam_boot_v)
hist(lam_boot_v, freq=F)
quantile(lam_boot_v, prob=c(0.025,0.975) )


# 6. Uncertainty via CI ------------------------------------------

gr_mod  <- lm( x1 ~ x, data = demo_df )
fec_mod <- glm( r0 ~ x, data = demo_df, family='poisson' )
sr_mod  <- glm( s1 ~ 1, data = demo_df, family='binomial' )
rec_mod <- lm( rcr ~ 1) # model of the mean for recruits

confint(gr_mod)

# lower parameters 
pars_1  <- list( surv_b0 = confint(sr_mod)[1],
               
                 grow_b0 = confint(gr_mod)[1,1],
                 grow_b1 = confint(gr_mod)[2,1],
                 grow_sd = summary(gr_mod)$sigma,
                 
                 fecu_b0 = confint(fec_mod)[1,1],
                 fecu_b1 = confint(fec_mod)[2,1],
                 
                 recr_sz = coef(rec_mod),
                 recr_sd = summary(rec_mod)$sigma,
                 L       = min(rcr,na.rm=T),
                 U       = max(c(demo_df$x,
                                 demo_df$x1),
                               na.rm=T),
                 mat_siz = 50
              )

pars_2  <- list( surv_b0 = confint(sr_mod)[2],
                 
                 grow_b0 = confint(gr_mod)[1,2],
                 grow_b1 = confint(gr_mod)[2,2],
                 grow_sd = summary(gr_mod)$sigma,
                 
                 fecu_b0 = confint(fec_mod)[1,2],
                 fecu_b1 = confint(fec_mod)[2,2],
                 
                 recr_sz = coef(rec_mod),
                 recr_sd = summary(rec_mod)$sigma,
                 L       = min(rcr,na.rm=T),
                 U       = max(c(demo_df$x,
                                 demo_df$x1),
                               na.rm=T),
                 mat_siz = 50
)


Re(eigen(kernel(pars_1)$k_yx)$value[1])
Re(eigen(kernel(pars_2)$k_yx)$value[1])


# 7. Uncertainty via Bayes ------------------------------------------

# set rstan options
rstan_options( auto_write = TRUE )
options( mc.cores = parallel::detectCores() )

# simulation parameters
sim_pars <- list(
  warmup = 1000, 
  iter = 4000, 
  thin = 2, 
  chains = 3
)

# growth model
gr_dat  <- list(
  n  = length(x),
  x  = demo_df$x,
  y  = demo_df$x1 
)

# fit growth model
gr_mod  <- stan(
  file = 'R/stan/gr_model.stan',
  data = gr_dat,
  pars = c('beta0','beta1','y_sd'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains
)

# growth model
fec_dat  <- list(
  n  = length(x),
  x  = demo_df$x,
  y  = demo_df$r0 
)

# fit growth model
fec_mod  <- stan(
  file = 'R/stan/pois_model.stan',
  data = fec_dat,
  pars = c('beta0','beta1'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains
)

# Surv model
sr_dat  <- list(
  n  = length(x),
  y  = demo_df$s1 
)

# fit growth model
sr_mod  <- stan(
  file = 'R/stan/bernoulli_model.stan',
  data = sr_dat,
  pars = c('beta0'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains
)

# posteriors
gr_p  <- rstan::extract(gr_mod) %>% as.data.frame
fec_p <- rstan::extract(fec_mod) %>% as.data.frame
sr_p  <- rstan::extract(sr_mod) %>% as.data.frame

# produce posterior estimates for lambda
post_lam <- function(ii){

  par_post <- list(  surv_b0 = sr_p[ii,'beta0'],
               
                     grow_b0 = gr_p[ii,'beta0'],
                     grow_b1 = gr_p[ii,'beta1'],
                     grow_sd = gr_p[ii,'y_sd'],
                     
                     fecu_b0 = fec_p[ii,'beta0'],
                     fecu_b1 = fec_p[ii,'beta1'],
                     
                     recr_sz = coef(rec_mod),
                     recr_sd = summary(rec_mod)$sigma,
                     L       = min(rcr,na.rm=T),
                     U       = max(c(demo_df$x,
                                     demo_df$x1),
                                   na.rm=T),
                     mat_siz = 50
              )
  
  Re(eigen(kernel(par_post)$k_yx)$value[1])
  
}

# sequence of numbers to "thin" the posterior
id_v <- round(seq(1,4500,length.out=1000))
id_v <- 1:4500

# vector of posterior estimates of lambda
lam_post_v <- sapply(id_v, post_lam)
quantile(lam_post_v, prob = c(0.025,0.975) )
