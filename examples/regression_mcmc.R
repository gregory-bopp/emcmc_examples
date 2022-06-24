library(tidyverse)
library(here)
library(purrr)
library(dplyr)
library(lobstr)        # For object sizes
library(microbenchmark)

set.seed(1)
##
## Simulate data
##
n <- 500
beta <- c(1, 100)
sigma_sq_y <- 0.5
X <- cbind(1, seq(0, 1, l = n))
y <- X %*% beta + rnorm(n, sd = sqrt(sigma_sq_y))
data <- list(X = X,
             y = y)
##
## Priors
##
priors <- list(
  beta = Normal$new(mean = c(0, 0),
                    sd = 100),
  sigma_sq_y = InvGamma$new(ig_shape = 0.01,
                            ig_scale = 1)
)

##
## Define likelihood
##
reg_ll_fn <- function(param,
                      data) {
  # param is a named list with elements
  # beta and sigma_sq_y
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

##
## Define and Run MCMC
##
set.seed(1)
# Setup MCMC
sampler <- MCMC$new(
  data = data,
  priors = priors,
  init = list(beta = beta,
              sigma_sq_y = sigma_sq_y),
  const = list(sigma_sq_y = 0.55),
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

