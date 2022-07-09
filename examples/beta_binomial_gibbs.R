library(tidyverse)
library(here)
library(purrr)
library(dplyr)
library(lobstr)        # For object sizes
library(microbenchmark)
library(zeallot)
library(emcmc)
# devtools::load_all(here('emcmc'))

set.seed(1)

# Simulate Data -----------------------------------------------------------
r_obs <- function() {
  n <- 500
  N <- 50
  shape1 <- 9
  shape2 <- 2
  theta <- rbeta(N, shape1 = shape1, shape2 = shape2)
  x <- rbinom(n = N, size = n, prob = theta)
  return(list(
    data = list(x = x,
                n = n,
                N = N),
    param = list(
      theta = theta,
      shape1 = shape1,
      shape2 = shape2
    )
  ))
}
c(data, gt) %<-% r_obs()



# Define Likelihood -------------------------------------------------------
# x | theta
ll_fn <- function(param,
                  data) {
  ll <-
    sum(dbinom(
      data$x,
      size = data$n,
      prob = param$theta,
      log = T
    ))

  if (is.na(ll))
    ll <- -Inf
  ll
  invisible(ll)
}

# theta|shape1, shape2
ll_theta <- function(param,
                     data) {
  ll <- sum(dbeta(
    x = param$theta,
    shape1 = param$shape1,
    shape2 = param$shape2,
    log = T
  ))
  if (is.na(ll))
    ll <- -Inf
  ll
  invisible(ll)
}

ll <-
  LogLik$new(cll_list = list(
    CompLogLik$new(param_names = c("theta"),
                   fn = ll_fn),
    CompLogLik$new(
      param_names = c("theta", "shape1", "shape2"),
      fn = ll_theta
    )
  ))



# Priors ------------------------------------------------------------------
priors <- list(shape1 = Normal$new(mean = 1,
                                   sd = 0.5),
               shape2 = Normal$new(mean = 1,
                                   sd = 0.5))


# Gibbs Updates -----------------------------------------------------------
theta_gibbs <- Proposal$new(
  update_type = "gibbs",
  r_fn = function(self) {
    value <- with(self$mcmc$cur, {
      value <- rbeta(data$N,
                     shape1 + self$mcmc$data$x,
                     shape2 + self$mcmc$data$n - data$x)
    })
    return(value)
  }
)


# Define and Run MCMC -----------------------------------------------------
# options(error=recover)
set.seed(1)
# Setup MCMC
sampler <- MCMC$new(
  data = data,
  priors = NA,
  init = list(
    theta = gt$theta,
    shape1 = gt$shape1,
    shape2 = gt$shape2
  ),
  proposals = list(theta = theta_gibbs),
  log_lik = ll,
  n_mcmc = 5000,
  thin_int = 10,
  show_progress = T
)

# Run MCMC
st <- Sys.time()
sampler$run_mcmc()
tot <- Sys.time() - st
print(tot)


# Posterior Summary -------------------------------------------------------
(ps <- sampler$post_summary())

# True vs posterior mean
plot(gt$theta, ps$theta$mean)

## Plot posterior samples
# Histogram
plts <- sampler$plot_post_histogram()
cowplot::plot_grid(plotlist = plts, nrow = 2)
# Traceplots
plts <- sampler$plot_post_trace()
cowplot::plot_grid(plotlist = plts, nrow = 2)
