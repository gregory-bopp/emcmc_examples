library(tidyverse)
library(here)
library(purrr)
library(dplyr)
library(zeallot)
library(microbenchmark)

########################################################################
## Model
########################################################################
#      y|beta, gamma, sigma_sq_y ~ N(X%*%beta + Z%*%gamma, sigma_sq_y)
#      beta|mu, sigma_sq_bg ~ N(mu, sigma_sq_bg)
#      gamma|mu, sigma_sq_bg ~ N(mu, sigma_sq_bg)
#      mu ~ N(0, 10)
#      sigma_sq_y ~ IG(a_y,b_y)
#      sigma_sq_bg ~ IG(a_bg, b_bg)
########################################################################

##
## Simulate data
##
r_obs <- function(){
  n <- 1000
  mu <- c(0,10)
  sigma_sq_bg <- 1/rgamma(1, shape = 1, rate = 0.1)
  beta <- rnorm(2, mean = mu, sd = sqrt(sigma_sq_bg))
  gamma <- rnorm(2, mean = mu, sd = sqrt(sigma_sq_bg))
  sigma_sq_y <- 1/rgamma(1, shape = 1, rate = 0.1)
  X <- cbind(1, seq(0, 1, l = n))
  Z <- cbind(rnorm(n), rnorm(n))
  y <- X %*% beta + Z %*% gamma + rnorm(n, sd = sqrt(sigma_sq_y))
  return(list(data = list(X = X,
               Z = Z,
               y = y),
              param = list(mu = mu,
                           sigma_sq_bg = sigma_sq_bg,
                           beta = beta,
                           gamma = gamma,
                           sigma_sq_y = sigma_sq_y
                           )
  )
  )
}
c(data, gt_param) %<-% r_obs()


##
## Priors
##
priors <- list(
  mu = Normal$new(mean = c(0),
                  sd = 10),
  sigma_sq_y = InvGamma$new(ig_shape = 1,
                            ig_scale = 0.1),
  sigma_sq_bg = InvGamma$new(ig_shape = 1,
                             ig_scale = 0.1)
)

##
## Define likelihood
##
# [y|beta, gamma, sigma_sq_y]
reg_ll_fn <- function(param,
                      data) {
  # param is a named list with elements
  # beta and sigma_sq_y
  #
  ll <- with(param, {
    ll <- sum(dnorm(
      data$y - data$X %*% beta - data$Z %*% gamma,
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
# [beta|mu, sigma_sq_bg] * [gamma|mu, sigma_sq_bg]
com_mean_ll_fn <- function(param,
                           data) {
  ll <- with(param, {
    sum(dnorm(
      beta,
      mean = mu,
      sd = sqrt(sigma_sq_bg),
      log = T
    )) +
      sum(dnorm(
        gamma,
        mean = mu,
        sd = sqrt(sigma_sq_bg),
        log = T
      ))
  })
  if (is.na(ll)) {
    ll <- -Inf
  }
  invisible(ll)
}

# Factorized log-likelihood components
reg_cll <- CompLogLik$new(param_names = c("beta","gamma","sigma_sq_y"),
                      fn = reg_ll_fn)
cm_cll <-  CompLogLik$new(param_names = c("beta","gamma","mu", "sigma_sq_bg"),
                      fn = com_mean_ll_fn)
# Aggregate Log-likelihood
ll <- LogLik$new(cll_list = list(reg_cll, cm_cll))

##
## Initial Values
##
init <- list(
  beta = c(0,0),
  gamma = c(0,0),
  sigma_sq_y = 1,
  sigma_sq_bg = 1,
  mu = c(0,0)
)

##
## Define and Run MCMC
##
set.seed(1)
# Setup MCMC
sampler <- MCMC$new(
  data = data,
  priors = priors,
  init = init,
  # const = list(sigma_sq_y = sigma_sq_y),
  log_lik = ll,
  n_mcmc = 10000,
  thin_int = 10,
  n_adapt = 1000,
  n_burnin = 1000
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
