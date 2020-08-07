##' Create forecasts based on a
##'
##' Take a grid search produced by \code{\link{grid_search}} and
##' sample \code{n_sample_pars} from the parameter grid uses based
##' on their probability. For each parameter pair chosen, run particle
##' filter with \code{num_particles}.
##'
##' @title Sample Grid Scan
##'
##' @param x Output of \code{\link{grid_search}}.
##'
##' @param ... Other arguments
##'
##' @param filter A \code{particle_filter} object to run (the same
##' (as the one used to produce \code{scan_results})
##'
##' @param n_sample_pars Number of parameters to be sampled. For a
##'   \code{mcstate_scan} these will be pairs of parameters sampled in
##'   proportion to their posterior probability. For a \code{mcstate_pmcmc_list}
##'   this will be full parameter sets sampled evenly from the chains.
##'   Integer. This will determine how many trajectories are returned.
##'
##' @param forecast_steps Number of time steps being forecast. Default = 0
##'
##' @param burn_in Amount of burn in to remove (if an \code{mcstate_pmcmc_list})
##'
##' @return \code{\link{list}}. First element (trajectories) is a 3
##'   dimensional array of trajectories (time, state, tranjectories). Second
##'   element (param_grid) is the parameters chosen when sampling from the
##'   \code{scan_results} grid and the third dimension (inputs) is a list of
##'   model inputs.
##'
##' @export
##' @importFrom utils tail
forecast <- function(x, filter, n_sample_pars = 10, forecast_steps = 0,
                     burn_in = 1) {
  browser()
  assert_integer(n_sample_pars)
  assert_integer(forecast_steps)

  # discard burn-in and sample
  par_names <- colnames(x$chains[[1]]$results)
  par_names <- par_names[!(par_names %in% c("log_prior",
                                            "log_likelihood",
                                            "log_posterior"))]
  chains <- combine_mcmc_chains(x, burn_in)
  param_grid <- chains[sample.int(n = nrow(chains),
                                  size = n_sample_pars,
                                  replace = FALSE), par_names]

  traces <- lapply(split_df_rows(param_grid), run_and_forecast,
                   filter, x$vars$index, forecast_steps)

  # combine and return
  res <- list("trajectories" = traces,
              "parameters" = param_grid)

  class(res) <- "mcstate_forecast"
  res

}

run_and_forecast <- function(pars, filter, forecast_steps) {
  filter$run(pars, save_history = TRUE)
  if (forecast_steps > 0) {
    forward_steps <- seq.int(0, forecast_steps)
    trajectories <- filter$predict(forward_steps, append = TRUE)
  } else {
    trajectories <- filter$history
  }
  trajectories
}
