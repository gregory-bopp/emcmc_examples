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

source(here("examples/sang_gelfand_2010/functions/block_gev_loc.R"))
source(here("examples/sang_gelfand_2010/functions/evd_functions.R"))

RFoptions(spConform = F)
set.seed(1)

# Data simulation function ------------------------------------------------
r_obs <- function(n_rep = 50,
                  x = seq(0, 10, len = 20),
                  gp_scale = 3,
                  gp_smooth = 1.5,
                  gev_loc_gp_scale = 4,
                  gev_loc_gp_var = 1,
                  gev_loc_gp_smooth = 1.5,
                  gev_scale = 1,
                  gev_shape = 0){
  n_loc <- length(x)
  # GEV location Covariance Function
  gev_loc_model <- RandomFields::RMmatern(
    var = 1,
    scale = gev_loc_gp_scale,
    nu = gev_loc_gp_smooth,
    notinvnu = T
  )
  gev_loc <- RandomFields::RFsimulate(n = 1,
                                     model = gev_loc_model,
                                     x = x) %>% c

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
  Y <- qgev(p = pnorm(Z),
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
      gev_loc_gp_var = gev_loc_gp_var,
      gev_loc_gp_scale = gev_loc_gp_scale,
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
  gev_scale = Normal$new(mean = 0, sd = 10),
  gev_shape = Normal$new(mean = 0, sd = 1),
  gev_loc_gp_scale = Normal$new(mean = 0, sd = 5),
  gev_loc_gp_var = Normal$new(mean = 0, sd = 5)
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
      Z <- qnorm(pgev(data$Y,
                      loc = gev_loc,
                      scale = gev_scale,
                      shape = gev_shape))

      ll <- RandomFields::RFlikelihood(model = model,
                                       x = data$x,
                                       data = Z)$loglikelihood
      # Jacobian Term
      ll <- ll + sum(dgev(data$Y,
                                    loc = gev_loc,
                                    scale = gev_scale,
                                    shape = gev_shape,
                                    log = T) -
                       dnorm(Z, sd = 1, log = T))
    }
  })
  if (is.na(ll)) {
    ll <- -Inf
  }
  invisible(ll)
}

gev_loc_ll <-  function(param, data) {
    ll <- with(param, {
      if (gev_loc_gp_scale <= 0 || gev_loc_gp_var <= 0) {
        ll <- -Inf
      }
      else{
        gev_loc_model <- RandomFields::RMmatern(var = gev_loc_gp_var,
                                        scale = gev_loc_gp_scale,
                                        nu = 1.5,
                                        notinvnu = T)

        ll <- RandomFields::RFlikelihood(model = gev_loc_model,
                                         x = data$x,
                                         data = gev_loc)$loglikelihood
      }
    })
    if (is.na(ll)) {
      ll <- -Inf
    }
    invisible(ll)
  }


# Component log-likelihoods
obs_cll <-
  CompLogLik$new(
    param_names = c("gev_loc", "gev_scale", "gev_shape",
                    "gp_scale", "gp_smooth"),
    fn = obs_ll
  )
gev_loc_cll <-
  CompLogLik$new(
    param_names = c("gev_loc_gp_scale", "gev_loc_gp_var"),
    fn = gev_loc_ll
  )
# Aggregate Log-likelihood
loglik <- LogLik$new(cll_list = list(obs_cll, gev_loc_cll))


# Initial Values
init <- list(gev_loc = gt$gev_loc,
             gev_loc_gp_var = gt$gev_loc_gp_var,
             gev_loc_gp_scale = gt$gev_loc_gp_scale,
             gev_scale = 5,
             gev_shape = 0.5,
             gp_scale = 2,
             gp_smooth = gt$gp_smooth)
init <- gt

proposals <- list(gev_loc = BlockGEVLoc$new(prop_var = 0.001,
                                            blocks = factor(rep(1, length(gt$gev_loc)))))

# Setup and Run MCMC ------------------------------------------------------
set.seed(1)
# Instantiate MCMC sampler
sampler <- MCMC$new(
  data = data,
  priors = priors,
  init = init,
  const = list(gp_smooth = gt$gp_smooth),
  proposals = proposals,
  log_lik = loglik,
  n_mcmc = 10000,
  thin_int = 10,
  n_adapt = 10000,
  n_burnin = 2000
)

# Run MCMC
sampler$run_mcmc()


# Posterior Summary -------------------------------------------------------
(ps <- sampler$post_summary())

# Histogram
plts <- sampler$plot_post_histogram()
cowplot::plot_grid(plotlist = plts, nrow = 3)
cowplot::plot_grid(plotlist = plts[names(plts) != "gev_loc"], nrow = 2)
plts[names(plts) == "gev_loc"]


# Traceplots
plts <- sampler$plot_post_trace()
cowplot::plot_grid(plotlist = plts[names(plts) != "gev_loc"], nrow = 2)
plts[names(plts) == "gev_loc"]


# Plot GEV location ------------------------------------------------------
ps_gev_loc <- ps$gev_loc %>%
  bind_cols(tibble(x = data$x))
gt_gev_loc <- data.frame(gt = gt$gev_loc, x = data$x)
ggplot(ps_gev_loc) +
  geom_ribbon(aes(x = x, ymin = `2.5%`, ymax = `97.5%`), alpha = 0.5) +
  geom_line(aes(x = x, y = mean)) +
  geom_line(aes(x = x, y = gt), col = 2,  data = gt_gev_loc)

