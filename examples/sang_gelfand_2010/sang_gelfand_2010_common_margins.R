#
# Description: MCMC implementation of the Sang and Gelfand (2010) model
#
#     Sang, H., & Gelfand, A. E. (2010). Continuous spatial process models
#     for spatial extreme values. Journal of agricultural, biological,
#     and environmental statistics, 15(1), 49-65.
#
library(tidyverse)
library(here)
library(purrr)
library(dplyr)
library(zeallot)
library(microbenchmark)
library(RandomFields)
library(extRemes)
library(emcmc)

rm(list = ls())

RFoptions(spConform = F)
set.seed(1)

# Data simulation function ------------------------------------------------
r_obs <- function(n_rep = 100,
                  x = seq(0, 10, len = 20),
                  gp_scale = 3,
                  gp_smooth = 1.5,
                  gev_loc = 0,
                  gev_scale = 1,
                  gev_shape = 0){
  n_loc <- length(x)
  # Correlation function
  model <- RandomFields::RMmatern(
    var = 1,
    scale = gp_scale,
    nu = gp_smooth,
    notinvnu = T
  )
  Z <- RandomFields::RFsimulate(n = n_rep,
                                model = model,
                                x = x)
  Y <- qevd(
    p = pnorm(Z),
    loc = gev_loc,
    scale = gev_scale,
    shape = gev_shape)

  return(list(
    data = list(x = x,
                Y = Y,
                Z = Z),
    param = list(
      gp_scale = gp_scale,
      gp_smooth = gp_smooth,
      gev_loc = gev_loc,
      gev_scale = gev_scale,
      gev_shape = gev_shape
    )
  ))
}

# Simulate data
c(data, gt) %<-% r_obs()



# Priors ------------------------------------------------------------------
priors <- list(
  gp_scale = Normal$new(mean = 0, sd = 5),
  gp_smooth = Normal$new(mean = 0, sd = 0.75),
  gev_loc = Normal$new(mean = 0, sd = 10),
  gev_scale = Normal$new(mean = 0, sd = 10),
  gev_shape = Normal$new(mean = 0, sd = 1)
)


# Likelihoods -------------------------------------------------------------
# Log-likelihood function for GP-GEV
obs_ll <- function(param, data) {
  ll <- with(param, {
    if (gp_scale <= 0 || gp_smooth <= 0 || gev_scale <= 0) {
      ll <- -Inf
    }
    else{
      model <- RandomFields::RMmatern(var = 1,
                                      scale = gp_scale,
                                      nu = gp_smooth,
                                      notinvnu = T)
      Z <- qnorm(extRemes::pevd(data$Y,
                                loc = gev_loc,
                                scale = gev_scale,
                                shape = gev_shape,
                                type = "GEV",
                                lower.tail = T,
                                log.p = F))
      ll <- RandomFields::RFlikelihood(model = model,
                                       x = data$x,
                                       data = Z)$loglikelihood
      # Jacobian Term
      ll <- ll + sum(extRemes::devd(data$Y,
                                    loc = gev_loc,
                                    scale = gev_scale,
                                    shape = gev_shape,
                                    type = "GEV",
                                    log = T) -
                       dnorm(Z, sd = 1, log = T))
    }
  })
  if (is.na(ll)) {
    ll <- -Inf
  }
  invisible(ll)
}

# Component log-likelihoods
obs_cll <- CompLogLik$new(param_names = c("gev_loc", "gev_scale", "gev_shape",
                                           "gp_scale", "gp_smooth"),
                          fn = obs_ll)
# Aggregate Log-likelihood
loglik <- LogLik$new(cll_list = list(obs_cll))



# Initial Values
init <- list(gev_loc = 3,
             gev_scale = 5,
             gev_shape = 0.5,
             gp_scale = 2,
             gp_smooth = gt$gp_smooth)
init <- gt


# Setup and Run MCMC ------------------------------------------------------
set.seed(1)
# Instantiate MCMC sampler
sampler <- MCMC$new(
  data = data,
  priors = priors,
  init = init,
  const = list(gp_smooth = gt$gp_smooth),
  log_lik = loglik,
  n_mcmc = 10000,
  thin_int = 10,
  n_adapt = 5000,
  n_burnin = 1000
)

# Run MCMC
sampler$run_mcmc()



# Posterior Summary -------------------------------------------------------
(ps <- sampler$post_summary())

# Histogram
plts <- sampler$plot_post_histogram()
cowplot::plot_grid(plotlist = plts, nrow = 3)

# Traceplots
plts <- sampler$plot_post_trace()
cowplot::plot_grid(plotlist = plts, nrow = 3)

