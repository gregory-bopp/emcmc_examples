require(zeallot)
require(magrittr)

BlockGEVLoc <- R6::R6Class(
  "BlockGEVLoc",
  inherit = Proposal,
  public =
    list(name = "BlockGEVLoc"),
  private =
    list(
      .adapt_prop_var = TRUE,
      .is_asymmetric = FALSE,
      .prop_var = NA
    )
)

BlockGEVLoc$set('public',
               'initialize',
               function(prop_var,
                        adapt_prop_var,
                        blocks) {
                 # Proposal variance
                 if (!missing(prop_var))
                   private$.prop_var <- prop_var
                 if (!missing(adapt_prop_var))
                   private$.adapt_prop_var <- adapt_prop_var
                 # Set blocks
                 self$set_blocks(blocks)
                 invisible(self)
               })

BlockGEVLoc$set('public',
               'r_fn',
               function(self, cur) {
                 c(gev_loc_gp_var, gev_loc_gp_scale) %<-% self$mcmc$cur[c("gev_loc_gp_var", "gev_loc_gp_scale")]
                 gev_loc_model <-
                   RandomFields::RMmatern(
                     var = gev_loc_gp_var,
                     scale = gev_loc_gp_scale,
                     nu = 1.5,
                     notinvnu = T
                   )
                 Cmat <-
                   RandomFields::RFcovmatrix(gev_loc_model, self$mcmc$data$x)
                 Lt <- chol(Cmat)
                 white <- cur %*% solve(Lt)
                 white_shifted <-
                   white + rnorm(length(cur),
                                 mean = 0,
                                 sd = sqrt(self$block_info$prop_var))
                 prop <- white_shifted %*% Lt %>% c
                 return(prop)
               })
