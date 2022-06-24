library(tidyverse)
library(here)
library(purrr)
library(dplyr)
library(zeallot)
library(microbenchmark)
library(RandomFields)
library(mvtnorm)
library(emcmc)

# RandomFields Options to return matrix instead of sp object from RFsimulate
RFoptions(spConform = F)

# Data simulation function ------------------------------------------------
r_obs <- function(n_rep = 50,
                  n_loc = 100,
                  gp_var = 1,
                  gp_scale = 3,
                  beta = c(1, 3),
                  gp_smooth = 1.5) {
  # Grid of spatial locations
  x = seq(0, 10, l = n_loc)
  # Distances between grid points
  D <- dist(x)
  # Simulate some known covariates
  Z <- cbind(seq(0, 10, l = n_loc),
             rnorm(n_loc))
  # Define the Gaussian Process Covariance model to use
  model <- RandomFields::RMmatern(
    var = gp_var,
    scale = gp_scale,
    nu = gp_smooth,
    notinvnu = TRUE
  )
  # Simulate the GP Errors
  eps <- map_dfc(1:n_rep, ~ {
    RandomFields::RFsimulate(n = 1,
                             model = model,
                             distances = D,
                             dim=c(n_loc))
  }) %>% as.matrix
  # GP Regression
  Y <- c(Z %*% beta) + eps

  # Return data and true parameters
  return(list(
    data = list(Y = Y, D = D, Z = Z),
    param = list(
      gp_var = gp_var,
      gp_scale = gp_scale,
      gp_smooth = gp_smooth,
      beta = beta
    )
  ))
}
# Simulate data
c(data, gt) %<-% r_obs()

# Define likelihood Object ------------------------------------------------
# Log-likelihood function for GP Regression
gp_ll <- function(param, data) {
  ll <- with(param, {
    n_loc <- nrow(data$Y)
    if (gp_var <= 0 || gp_scale <= 0 || gp_smooth <= 0) {
      ll <- -Inf
    }
    else{
      # Define covariance model
      model <- RandomFields::RMmatern(
        var = gp_var,
        scale = gp_scale,
        nu = gp_smooth,
        notinvnu = TRUE
      )
      # Calculate covariance matrix
      K <-  RandomFields::RFcovmatrix(model = model,
                                        distances = data$D,
                                        dim = n_loc)
      # Take Get GP Errors
      eps <-  data$Y - c(data$Z %*% beta)
      # Calculate log-likelihood
      ll <- sum(dmvnorm(t(eps), mean = rep(0, n_loc), sigma = K, log = T))
    }
    ll
  })
  if (is.na(ll)) {
    ll <- -Inf
  }
  invisible(ll)
}
# A component log-likelihood object (i.e. log of factorized likelihood components)
gp_cll <-
  CompLogLik$new(param_names = c("gp_var", "gp_scale", "gp_smooth", "beta"),
                 fn = gp_ll)
# Aggregate Log-likelihood
ll <- LogLik$new(cll_list = list(gp_cll))


# Priors ------------------------------------------------------------------
priors <- list(
  gp_var =  InvGamma$new(ig_shape = 1, ig_scale = 0.1),
  gp_scale = Normal$new(mean = 0, sd = 5),       # Half-normal (constrained below at zero in ll)
  gp_smooth = Normal$new(mean = 0, sd = 0.75),   # Half-normal (constrained below at zero in ll)
  beta = Normal$new(mean = 0, sd = 10)
)

##
## Initial Values
##
init <- list(gp_var = 5,
             gp_scale = 1,
             # gp_smooth = 0.5,
             beta = c(2, 5))


##
## Define and Run MCMC
##
options(error = NULL)
set.seed(1)
# Setup MCMC
sampler <- MCMC$new(
  data = data,
  priors = priors,
  init = init,
  const = list(gp_smooth = gt$gp_smooth),  # Constants (not sampled)
  log_lik = ll,
  n_mcmc = 5000,
  thin_int = 10,
  n_adapt = 3000,
  n_burnin = 500
)

# Run MCMC
start_t <- Sys.time()
sampler$run_mcmc()
tot_t <- Sys.time() - start_t
print(tot_t)

##
## Posterior Summary
##
sampler$post_summary()

## Plot posterior samples
# Histogram
plts <- sampler$plot_post_histogram()
cowplot::plot_grid(plotlist = plts, nrow = 2)
# Traceplots
plts <- sampler$plot_post_trace()
cowplot::plot_grid(plotlist = plts, nrow = 2)
