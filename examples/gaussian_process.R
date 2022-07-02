library(tidyverse)
library(here)
library(purrr)
library(dplyr)
library(zeallot)
library(microbenchmark)
library(RandomFields)
library(emcmc)

RFoptions(spConform = F)


# Data simulation function ------------------------------------------------
r_obs <- function(n_rep = 50,
                  x = seq(0, 10, len = 100),
                  gp_var = 1,
                  gp_scale = 3,
                  gp_smooth = 1.5) {
  model <- RandomFields::RMmatern(
    var = gp_var,
    scale = gp_scale,
    nu = gp_smooth,
    notinvnu = TRUE
  )
  Y <- map_dfc(1:n_rep, ~ {
    RandomFields::RFsimulate(n = 1,
                             model = model,
                             x = x)
  }) %>% as.matrix
  return(list(
    data = list(Y = Y, x = x),
    param = list(
      gp_var = gp_var,
      gp_scale = gp_scale,
      gp_smooth = gp_smooth
    )
  ))
}
# Simulate data
c(data, gt) %<-% r_obs()


# Log-likelihood function for GP ------------------------------------------
gp_ll <- function(param, data) {
  ll <- with(param, {
    if (gp_var <= 0 || gp_scale <= 0 || gp_smooth <= 0) {
      ll <- -Inf
    }
    else{
      model <- RandomFields::RMmatern(
        var = gp_var,
        scale = gp_scale,
        nu = gp_smooth,
        notinvnu = TRUE
      )
      ll <-  RandomFields::RFlikelihood(model = model,
                                        x = data$x,
                                        data = data$Y)$loglikelihood
    }
    ll
  })
  if (is.na(ll)) {
    ll <- -Inf
  }
  invisible(ll)
}


# Priors ------------------------------------------------------------------
priors <- list(
  gp_var =  Normal$new(mean = 0, sd = 10),
  gp_scale = Normal$new(mean = 0, sd = 5),
  gp_smooth = Normal$new(mean = 0, sd = 0.75)
)

gp_cll <-
  CompLogLik$new(param_names = c("gp_var", "gp_scale", "gp_smooth"),
                 fn = gp_ll)
# Aggregate Log-likelihood
ll <- LogLik$new(cll_list = list(gp_cll))



##
## Initial Values
##
init <- list(gp_var = 5,
             gp_scale = 1,
             gp_smooth = 0.5)


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
  const = list(gp_smooth = gt$gp_smooth),
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
