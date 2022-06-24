library(tidyverse)
library(here)
library(purrr)
library(dplyr)
library(microbenchmark)

##
## Base proposal to be copied for beta and sigma_sq
##
# Use self$prop_var to reference adaptively tuned prop_var
asym_prop <- Proposal$new(
  r_fn = function(self, cur) {
    return(rnorm(length(cur),
                 mean = cur + 0.1,
                 sd = 0.1))
  },
  d_fn = function(self, prop, cur) {
    return(dnorm(prop,
                 cur,
                 sd = 0.1,
                 log = T))
  },
  is_asymmetric = T,
  adapt_prop_var = T,
  prop_var = 0.1
)


##
## Simulate data
##
n <- 500
beta <- c(1, 0.5)
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

init <- list(beta = beta,
     sigma_sq_y = sigma_sq_y)



proposals <- map(names(init), ~{
  asym_prop$clone(deep = T)
  }) %>%
  `names<-`(names(init))

##
## Define and Run MCMC
##
set.seed(1)
# Setup MCMC
sampler <- MCMC$new(
  data = data,
  priors = priors,
  init = init,
  proposals = proposals,
  # const = list(sigma_sq_y = 0.55),
  log_lik = ll,
  n_mcmc = 10000,
  thin_int = 10
)
# Run MCMC
sampler$run_mcmc()


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


