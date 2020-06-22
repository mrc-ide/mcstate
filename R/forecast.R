##' Take a grid search produced by \code{\link{grid_search}} and
##' sample \code{n_sample_pairs} from the parameter grid uses based
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
##' @param n_sample_pairs Number of parameter pairs to be sampled. This will
##'   determine how many trajectories are returned. Integer. Default = 10. This
##'   will determine how many trajectories are returned.
##'
##' @param n_particles Number of particles. Positive Integer. Default = 100
##'
##' @param forecast_steps Number of time steps being forecast. Default = 0
##'
##' @return \code{\link{list}}. First element (trajectories) is a 3
##'   dimensional array of trajectories (time, state, tranjectories). Second
##'   element (param_grid) is the parameters chosen when sampling from the
##'   \code{scan_results} grid and the third dimension (inputs) is a list of
##'   model inputs.
##'
##' @export
##' @importFrom utils tail
forecast <- function(x, ...,
                     filter,
                     n_sample_pairs = 10,
                     n_particles = 100,
                     forecast_steps = 0) {
  UseMethod("forecast", x)
}

##' Forecast for grid search
##' @rdname forecast
##' @method forecast mcstate_scan
forecast.mcstate_scan <- function(x, ...,
                                  filter,
                                  n_sample_pairs = 10,
                                  n_particles = 100,
                                  forecast_steps = 0) {

  # checks on args
  assert_integer(n_sample_pairs)
  assert_integer(n_particles)

  # sample proportional to probability
  sample_idx <- sample(nrow(x$vars$expanded),
                       size = n_sample_pairs,
                       replace = TRUE,
                       prob = x$renorm_mat_ll)
  pairs <- x$vars$expanded[sample_idx, ]

  traces <- lapply(split_df_rows(pairs), run_and_forecast,
                   filter, x$vars$index, n_particles,
                   forecast_steps)

  # combine and return
  res <- list("trajectories" = traces,
              "parameters" = pairs)

  class(res) <- "mcstate_forecast"
  return(res)

}

run_and_forecast <- function(model_params, filter, index, n_particles,
                             forecast_steps) {
  filter$run2(n_particles, save_history = TRUE, index, model_params)
  if (forecast_steps > 0) {
    forward_steps <- seq.int(0, forecast_steps)
    trajectories <- filter$predict(forward_steps, append = TRUE)
  } else {
    trajectories <- filter$history
  }
  trajectories
}
