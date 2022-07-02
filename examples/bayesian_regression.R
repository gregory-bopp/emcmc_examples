library(tidyverse)
library(here)
library(purrr)
library(dplyr)
library(lobstr)        # For object sizes
library(microbenchmark)
library(zeallot)
library(emcmc)

set.seed(1)

# Simulate Data -----------------------------------------------------------
r_obs <- function() {
  n <- 500
  beta <- c(1, 100)
  sigma_sq_y <- 0.5
  X <- cbind(1, seq(0, 1, l = n))
  y <- X %*% beta + rnorm(n, sd = sqrt(sigma_sq_y))
  return(list(
    data = list(X = X,
                 y = y),
    param = list(beta = beta,
                 sigma_sq_y = sigma_sq_y)
  ))
}
c(data, gt) %<-% r_obs()

# Define Priors -----------------------------------------------------------
priors <- list(
  beta = Normal$new(mean = c(0, 0),
                    sd = 100),
  sigma_sq_y = InvGamma$new(shape = 0.01,
                            scale = 1)
)


# Define Likelihood -------------------------------------------------------
reg_ll_fn <- function(param,
                      data) {
  ll <- with(param, {
    ll <- sum(dnorm(
      data$y - data$X %*% beta,
      mean = 0,
      sd = sqrt(sigma_sq_y),
      log = T
    ))
    if (is.na(ll))
      ll <- -Inf
    ll
  })
  invisible(ll)
}
cll <- CompLogLik$new(param_names = c("beta", "sigma_sq_y"),
                      fn = reg_ll_fn)
ll <- LogLik$new(cll_list = list(cll))


# Define and Run MCMC -----------------------------------------------------
set.seed(1)
# Setup MCMC
sampler <- MCMC$new(
  data = data,
  priors = priors,
  init = list(beta = gt$beta,
              sigma_sq_y = gt$sigma_sq_y),
  # const = list(sigma_sq_y = gt$sigma_sq_y),
  log_lik = ll,
  n_mcmc = 5000,
  thin_int = 10,
  # cache_freq = 1000,
  # combine_cached = T,
  show_progress = T
)

# Run MCMC
st <- Sys.time()
sampler$run_mcmc()
tot <- Sys.time() - st
print(tot)


# Posterior Summary -------------------------------------------------------
sampler$post_summary()
## Plot posterior samples
# Histogram
plts <- sampler$plot_post_histogram()
cowplot::plot_grid(plotlist = plts, nrow = 2)
# Traceplots
plts <- sampler$plot_post_trace()
cowplot::plot_grid(plotlist = plts, nrow = 2)
