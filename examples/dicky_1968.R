library(tidyverse)
library(here)
library(purrr)
library(dplyr)
library(microbenchmark)
library(zeallot)
library(emcmc)



# Simulate Data -----------------------------------------------------------
r_obs <- function(){
  n <- 20
  beta <- c(1, 0.5)
  sigma_sq_theta <- 0.1
  sigma_sq_y <- 0.3
  X <- cbind(1, seq(0, 1, l = n))
  theta <- X %*% beta + rnorm(n, sd = sqrt(sigma_sq_theta))
  y <- rnorm(n, mean = theta, sd = sqrt(sigma_sq_y))
  return(list(data = list(X = X,
                           y = y),
              param = list(theta = theta,
                           beta = beta,
                           sigma_sq_theta = sigma_sq_theta,
                           sigma_sq_y = sigma_sq_y)
  ))
}
c(data, gt) %<-% r_obs()


# Define priors -----------------------------------------------------------
priors <- list(
  beta = Normal$new(mean = c(0, 0),
                    sd = 100),
  sigma_sq_y = InvGamma$new(shape = 1,
                            scale = 1),
  sigma_sq_theta = InvGamma$new(shape = 1,
                                scale = 1)
)

##
## Define likelihood
##
theta_ll <- function(param,
                      data) {
  ll <- with(param, {
    ll <- sum(dnorm(
      theta - data$X %*% beta,
      mean = 0,
      sd = sqrt(sigma_sq_theta),
      log = T
    ))
    if (is.na(ll))
      ll <- -Inf
    ll
  })
  invisible(ll)
}
obs_ll <- function(param,
                   data){
  ll <- with(param,{
    ll <- sum(dnorm(data$y, theta, sd = sqrt(sigma_sq_y),
                    log = T))
    if (is.na(ll))
      ll <- -Inf
    ll
  })
  invisible(ll)
}
# Component log-likelihood objects
cll_theta <- CompLogLik$new(param_names = c("beta", "theta", "sigma_sq_theta"),
                      fn = theta_ll)
cll_obs <- CompLogLik$new(param_names = c("theta", "sigma_sq_y"),
                          fn = obs_ll)
# Combine component log-likelihoods into log-likelihood
ll <- LogLik$new(cll_list = list(cll_theta,
                                 cll_obs))

##
## Define and Run MCMC --------------------------------------------------------
##
set.seed(1)
# Setup MCMC
sampler <- MCMC$new(
  data = data,
  priors = priors,
  init = gt,
  # const = list(sigma_sq_y = 0.55),
  log_lik = ll,
  n_mcmc = 10000,
  thin_int = 10
)

# Run MCMC
sampler$run_mcmc()
# Add mcmc
sampler$add_n_mcmc(1000)
sampler$run_mcmc()


# Posterior Summary -------------------------------------------------------
sampler$post_summary()
## Plot posterior samples
# Histogram
plts <- sampler$plot_post_histogram()
cowplot::plot_grid(plotlist = plts[names(plts) != "theta"], nrow = 3)
plts[names(plts) == "theta"]

# Traceplots
plts <- sampler$plot_post_trace()
cowplot::plot_grid(plotlist = plts[names(plts) != "theta"], nrow = 3)
plts[names(plts) == "theta"]


