##' Take a grid search produced by \code{\link{grid_search}} and
##' sample \code{n_sample_pairs} from the parameter grid uses based
##' on their probability. For each parameter pair chosen, run particle
##' filter with \code{num_particles}.
##'
##' @title Sample Grid Scan
##'
##' @param scan_results Output of \code{\link{grid_search}}.
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
##' @import purrr
##' @importFrom utils tail
sample_grid_scan <- function(scan_results,
                             filter,
                             n_sample_pairs = 10,
                             n_particles = 100,
                             forecast_steps = 0) {

  # checks on args
  assert_is(scan_results, "mcstate_scan")
  assert_integer(n_sample_pairs)
  assert_integer(n_particles)

  # sample proportional to probability
  sample_idx <- sample(nrow(scan_results$vars$expanded),
                       size = n_sample_pairs,
                       replace = TRUE,
                       prob = scan_results$renorm_mat_ll)
  pairs <- scan_results$vars$expanded[sample_idx, ]

  traces <- purrr::map(.x = purrr::transpose(pairs), .f = run_and_forecast,
                       filter, scan_results$vars$index, n_particles,
                       forecast_steps)

  # If the start point was sampled, trajectories will have different
  # lengths and need to be filled with NAs
  if (length(scan_results$vars$index$step_start) > 1) {
    traces <- traces_to_trajectories(traces)
  }

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

##' collapse into an array of trajectories
##' the trajectories are different lengths in terms of dates
##' so we will fill the arrays with NAs where needed
##' @importFrom utils tail
##' @noRd
traces_to_trajectories <- function(traces) {

  num_rows <- unlist(lapply(traces, nrow))
  max_rows <- max(num_rows)
  seq_max <- seq_len(max_rows)
  max_date_names <- rownames(traces[[which.max(unlist(lapply(traces, nrow)))]])
  trajectories <- array(NA,
                        dim = c(max_rows, ncol(traces[[1]]), length(traces)),
                        dimnames = list(max_date_names, NULL, NULL))

  # fill the tail of the array slice
  # This is so that the end of the trajectories array is populated,
  # and the start is padded with NA if it's shorter than the max.
  for (i in seq_len(length(traces))) {
    trajectories[tail(seq_max, nrow(traces[[i]])), , i] <- traces[[i]]
  }

  trajectories
}
