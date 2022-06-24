library(tidyverse)
library(here)
library(purrr)
library(dplyr)
library(zeallot)
library(microbenchmark)
library(RandomFields)
library(extRemes)

RFoptions(spConform = F)

set.seed(1)
# Data simulation function ------------------------------------------------
r_obs <- function(n_rep = 50,
                  x = seq(0, 10, len = 20),
                  gp_var = 1,
                  gp_scale = 3,
                  gp_smooth = 1.5,
                  gev_scale = 1,
                  gev_shape = 0) {
  model <- RandomFields::RMmatern(
    var = gp_var,
    scale = gp_scale,
    nu = gp_smooth,
    notinvnu = TRUE
  )
  n_loc <- length(x)
  gev_loc <-  RandomFields::RFsimulate(n = 1,
                             model = model,
                             x = x)
  Y <- map_dfc(1:n_rep, ~ {
      revd(n_loc, gev_loc, scale = gev_scale, shape = gev_shape)
  }) %>%
    as.matrix
  return(list(
    data = list(Y = Y, x = x),
    param = list(
      gp_var = gp_var,
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
  gp_var =  Normal$new(mean = 0, sd = 10),
  gp_scale = Normal$new(mean = 0, sd = 5),
  gp_smooth = Normal$new(mean = 0, sd = 0.75),
  gev_loc = Normal$new(mean = 0, sd = 100),
  gev_scale = Normal$new(mean = 0, sd = 10),
  gev_shape = Normal$new(mean = 0, sd = 1)
)

# Likelihoods -------------------------------------------------------------
# Log-likelihood function for GP GEV Loc
gev_loc_ll <- function(param, data) {
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
                                        data = gev_loc)$loglikelihood
    }
    ll
  })
  if (is.na(ll)) {
    ll <- -Inf
  }
  invisible(ll)
}
# Observation Likelihoods
obs_ll <- function(param, data) {
  ll <- with(param, {
    extRemes::devd(data$Y)
    sum(apply(data$Y, 2, function(y)
      extRemes::devd(
        y,
        loc = gev_loc,
        scale = gev_scale,
        shape = gev_shape,
        log = T
      )))
  })
  if (is.na(ll)) {
    ll <- -Inf
  }
  invisible(ll)
}

# Create Componenet Likelihoods
gev_loc_cll <-
  CompLogLik$new(param_names = c("gev_loc", "gp_var", "gp_scale", "gp_smooth"),
                 fn = gev_loc_ll)
obs_cll <-
  CompLogLik$new(param_names = c("gev_loc", "gev_scale", "gev_shape"),
                 fn = obs_ll)
# Aggregate Log-likelihood
ll <- LogLik$new(cll_list = list(gev_loc_cll, obs_cll))


##
## Proposals
##
gev_loc_prop <- NormalRW$new(prop_var = 0.001,
                         blocks = factor(rep(1, length(data$x))))

##
## Initial Values
##
init <- list(gev_loc = rep(0, length(data$x)),
             gp_var = 5,
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
  const = list(gp_smooth = gt$gp_smooth,
               gev_scale = gt$gev_scale,
               gev_shape = gt$gev_shape),
  proposals = list(gev_loc = gev_loc_prop),
  log_lik = ll,
  n_mcmc = 10000,
  thin_int = 10,
  n_adapt = 2000,
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
(ps <- sampler$post_summary())

## Plot posterior samples
# Histogram
plts <- sampler$plot_post_histogram()
cowplot::plot_grid(plotlist = plts[names(plts) != "gev_loc"], nrow = 2)
plts[names(plts) == "gev_loc"]

# Traceplots
plts <- sampler$plot_post_trace()
cowplot::plot_grid(plotlist = plts[names(plts) != "gev_loc"], nrow = 3)
plts[names(plts) == "gev_loc"]



# Create Data for a Posterior Summary -------------------------------
df <- tibble(Y = c(data$Y),
             x = rep(data$x, ncol(data$Y)),
             rep = rep(1:ncol(data$Y), each = nrow(data$Y)))
par_df <- data.frame(x = data$x,
                     gt_gev_loc = gt$gev_loc,
                     mean_gev_loc = ps$gev_loc$mean,
                     lb = ps$gev_loc$`2.5%`,
                     ub = ps$gev_loc$`97.5%`)


# Plot data with inferred GEV location ------------------------------------
plist <- list(
  ggplot(df) +
    geom_point(aes(x = jitter(x,factor = 0.5), y = Y), size = 1) +
    labs(x = "x", y = "obs") +
    theme_bw(),
  ggplot(par_df) +
    geom_ribbon(aes(
      x = x, ymin = lb, ymax = ub
    )) +
    geom_line(aes(x = x, y = mean_gev_loc)) +
    geom_line(aes(x = x, y = gt_gev_loc), color = 'red') +
    labs(y = "GEV Location") +
    theme_bw()
)

cowplot::plot_grid(plotlist = plist, nrow = 2)

