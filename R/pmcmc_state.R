pmcmc_state <- R6::R6Class(
  "pmcmc_state",

  private = list(
    filter = NULL,
    pars = NULL,
    control = NULL,
    deterministic = NULL,

    nested = NULL,

    history_pars = NULL,
    history_probabilities = NULL,
    history_state = NULL,
    history_restart = NULL,
    history_trajectories = NULL,

    curr_step = NULL,
    curr_pars = NULL,
    curr_lprior = NULL,
    curr_llik = NULL,
    curr_lpost = NULL,
    curr_trajectories = NULL,
    curr_state = NULL,
    curr_restart = NULL,

    tick = NULL,
    update = NULL,
    adaptive = NULL,

    update_particle_history = function() {
      if (private$deterministic) {
        i <- 1L
      } else {
        i <- sample.int(private$filter$n_particles, 1)
      }

      if (private$control$save_trajectories) {
        private$curr_trajectories <- array_drop(private$filter$history(i), 2)
      }
      if (private$control$save_state) {
        private$curr_state <-
          array_drop(array_nth_dimension(private$filter$state(), 2, i), 2)
      }
      if (length(private$control$save_restart) > 0) {
        private$curr_restart <- array_drop(private$filter$restart_state(i), 2)
      }
    },

    update_mcmc_history = function(i) {
      private$history_pars$add(i, private$curr_pars)

      if (private$nested) {
        p <- rbind(private$curr_lprior, private$curr_llik, private$curr_lpost,
                   deparse.level = 0)
      } else {
        p <- c(private$curr_lprior, private$curr_llik, private$curr_lpost)
      }
      private$history_probabilities$add(i, p)

      control <- private$control
      i <- i - control$n_burnin - 1
      if (i >= 0 && i %% control$n_steps_every == 0) {
        j <- i / control$n_steps_every + 1
        if (!is.null(private$history_trajectories)) {
          private$history_trajectories$add(j, private$curr_trajectories)
        }
        if (!is.null(private$history_state)) {
          private$history_state$add(j, private$curr_state)
        }
        if (!is.null(private$history_restart)) {
          private$history_restart$add(j, private$curr_restart)
        }
      }
    },

    ## Computing the acceptance thresold, where u is a random uniform
    ## draw:
    ##
    ## => u < exp(prop_llik + prop_lprior - curr_lpost)
    ## => curr_lpost - prop_lprior + log(u) < prop_llik
    min_log_likelihood = function(prop_lprior, u) {
      if (!private$control$filter_early_exit) {
        return(-Inf)
      }
      if (length(prop_lprior) == length(u)) {
        private$curr_lpost - prop_lprior + log(u)
      } else {
        sum(private$curr_lpost - prop_lprior) + log(u)
      }
    },

    run_filter = function(p, min_log_likelihood = -Inf) {
      private$filter$run(private$pars$model(p),
                         private$control$save_trajectories,
                         private$control$save_restart,
                         min_log_likelihood)
    },

    update_simple = function() {
      browser()
      is_adaptive <- !is.null(private$adaptive)
      if (is_adaptive) {
        prop_pars <- private$adaptive$propose(private$curr_pars)
      } else {
        prop_pars <- private$pars$propose(private$curr_pars)
      }

      prop_lprior <- private$pars$prior(prop_pars)

      u <- runif(1)
      min_llik <- private$min_log_likelihood(prop_lprior, u)

      prop_llik <- private$run_filter(prop_pars, min_llik)
      prop_lpost <- prop_lprior + prop_llik

      accept <- u < exp(prop_lpost - private$curr_lpost)
      if (accept) {
        private$curr_pars <- prop_pars
        private$curr_lprior <- prop_lprior
        private$curr_llik <- prop_llik
        private$curr_lpost <- prop_lpost
        private$update_particle_history()
      }

      if (is_adaptive) {
        private$adaptive$update(private$curr_pars, accept)
      }
    },

    update_combined = function(type) {
      is_adaptive <- !is.null(private$adaptive)
      if (is_adaptive) {
        prop_pars <- private$adaptive$propose(private$curr_pars, type = type)
      } else {
        prop_pars <- private$pars$propose(private$curr_pars, type = type)
      }
      prop_lprior <- private$pars$prior(prop_pars)

      u <- runif(1)
      min_llik <- private$min_log_likelihood(prop_lprior, u)

      prop_llik <- private$run_filter(prop_pars, min_llik)
      prop_lpost <- prop_lprior + prop_llik

      accept <- u < exp(sum(prop_lpost - private$curr_lpost))
      if (accept) {
        private$curr_pars <- prop_pars
        private$curr_lprior <- prop_lprior
        private$curr_llik <- prop_llik
        private$curr_lpost <- prop_lpost
        private$update_particle_history()
      }

      if (is_adaptive) {
        private$adaptive$update(private$curr_pars, type = type, accept)
      }
    },

    update_fixed = function() {
      private$update_combined("fixed")
    },

    update_both = function() {
      private$update_combined("both")
    },

    update_varied = function() {
      type <- "varied"
      is_adaptive <- !is.null(private$adaptive)
      if (is_adaptive) {
        prop_pars <- private$adaptive$propose(private$curr_pars, type = type)
      } else {
        prop_pars <- private$pars$propose(private$curr_pars, type = type)
      }
      prop_lprior <- private$pars$prior(prop_pars)

      u <- runif(length(prop_lprior))
      min_llik <- private$min_log_likelihood(prop_lprior, u)

      prop_llik <- private$run_filter(prop_pars, min_llik)
      prop_lpost <- prop_lprior + prop_llik

      accept <- u < exp(prop_lpost - private$curr_lpost)
      if (any(accept)) {
        private$curr_pars[, accept] <- prop_pars[, accept]
        private$curr_lprior[accept] <- prop_lprior[accept]
        private$curr_llik[accept] <- prop_llik[accept]
        private$curr_lpost[accept] <- prop_lpost[accept]
        private$update_particle_history()
      }

      if (is_adaptive) {
        private$adaptive$update(private$curr_pars, type = type, accept)
      }
    }
  ),

  public = list(
     initialize = function(pars, initial, filter, control) {
      private$filter <- filter
      private$pars <- pars
      private$control <- control
      private$nested <- inherits(pars, "pmcmc_parameters_nested")
      private$deterministic <- inherits(filter, "particle_deterministic")

      if (filter$has_multiple_parameters && !filter$has_multiple_data) {
        stop(paste("Can't use a filter with multiple parameter sets but not",
                   "multiple data"))
      }

      if (private$nested != filter$has_multiple_data) {
        stop("'pars' and 'filter' disagree on nestedness")
      }

      if (!is.null(control$adaptive_proposal)) {
        if (!private$deterministic) {
          stop("Adaptive proposal only allowed in deterministic models")
        }
        if (private$nested) {
          private$adaptive <- adaptive_proposal_nested$new(
            pars, control$adaptive_proposal)
        } else {
          private$adaptive <- adaptive_proposal$new(
            pars, control$adaptive_proposal)
        }
      }

      private$tick <- pmcmc_progress(control$n_steps, control$progress,
                                     control$progress_simple)

      private$curr_step <- 0L
      private$curr_pars <- initial
      private$curr_lprior <- private$pars$prior(private$curr_pars)
      private$curr_llik <- private$run_filter(private$curr_pars)
      private$curr_lpost <- private$curr_lprior + private$curr_llik
      private$update_particle_history()

      n_steps <- control$n_steps
      n_history <- control$n_steps_retain

      private$history_pars <- history_collector(n_steps)
      private$history_probabilities <- history_collector(n_history)
      if (control$save_trajectories) {
        private$history_trajectories <- history_collector(n_history)
      }
      if (control$save_state) {
        private$history_state <- history_collector(n_history)
      }
      if (length(control$save_restart) > 0) {
        private$history_restart <- history_collector(n_history)
      }

      if (!private$nested) {
        update <- update_single(private$update_simple)
      } else if (length(pars$names("fixed")) == 0) {
        update <- update_single(private$update_varied)
      } else if (length(pars$names("varied")) == 0) {
        update <- update_single(private$update_fixed)
      } else if (private$control$nested_update_both) {
        update <- update_single(private$update_both)
      } else {
        update <- update_alternate(private$update_fixed,
                                   private$update_varied,
                                   private$control$nested_step_ratio)
      }
      private$update <- update
    },

    run = function() {
      control <- private$control
      ## TODO: simplify, then look at simplifying the rest
      to <- min(private$curr_step + control$n_steps, control$n_steps)
      steps <- seq(from = private$curr_step + 1L,
                   length.out = to - private$curr_step)
      rerun <- make_rerun(control$rerun_every, control$rerun_random)

      for (i in steps) {
        private$tick()

        if (rerun(i)) {
          private$curr_llik <- private$run_filter(private$curr_pars)
          private$curr_lpost <- private$curr_lprior + private$curr_llik
          private$update_particle_history()
        }

        private$update(i)
        private$update_mcmc_history(i)
      }

      private$curr_step <- to

      list(step = to, finished = to == control$n_steps)
    },

    finish = function() {
      nms_probabilities <- c("log_prior", "log_likelihood", "log_posterior")
      if (private$nested) {
        idx_pars <- c(3, 1, 2)
        idx_state <- c(1, 2, 4, 3)
        dimnames_pars <- c(list(NULL), dimnames(private$curr_pars))
        dimnames_probabilities <- list(NULL, nms_probabilities,
                                       private$pars$populations())
      } else {
        idx_pars <- c(2, 1)
        idx_state <- c(1, 3, 2)
        dimnames_pars <- list(NULL, names(private$curr_pars))
        dimnames_probabilities <- list(NULL, nms_probabilities)
      }

      ## sample x par | sample x par x pop
      pars <- array_from_list(
        private$history_pars$get(), idx_pars)
      dimnames(pars) <- dimnames_pars

      probabilities <- array_from_list(
        private$history_probabilities$get(), idx_pars)
      dimnames(probabilities) <- dimnames_probabilities

      predict <- state <- restart <- trajectories <- NULL

      if (private$control$save_state || private$control$save_trajectories) {
        ## TODO: tidy up private access here; check what uses this?
        ##
        ## Do we *definitely* need step and rate here?
        data <- private$filter$inputs()$data
        is_continuous <- inherits(data, "particle_filter_data_continuous")
        if (is_continuous) {
          step <- NULL
          rate <- NULL
          time <- last(data$time_end)
        } else {
          step <- last(data$step_end)
          rate <- attr(data, "rate", exact = TRUE)
          time <- step * rate
        }

        predict <- list(
          is_continuous = is_continuous,
          transform = r6_private(private$pars)$transform,
          index = r6_private(private$filter)$last_history$index,
          step = step,
          rate = rate,
          time = time,
          filter = private$filter$inputs())
      } else {
        predict <- NULL
      }

      if (private$control$save_state) {
        ## state x sample | state x pop x sample
        state <- array_from_list(private$history_state$get())
      }

      if (length(private$control$save_restart) > 0) {
        ## [state x sample x time] (from [state x time] x sample)
        ## [state x pop x sample x time] (from [state x pop x time] x sample)
        restart_state <-
          array_from_list(private$history_restart$get(), idx_state)
        restart <- list(time = private$control$save_restart,
                        state = restart_state)
      }

      if (private$control$save_trajectories) {
        ## [state x sample x time] (from [state x time] x sample)
        ## [state x pop x sample x time] (from [state x pop x time] x sample)
        trajectories_state <-
          array_from_list(private$history_trajectories$get(), idx_state)
        rownames(trajectories_state) <- names(predict$index)
        if (private$nested) {
          colnames(trajectories_state) <- private$pars$populations()
        }
        ## This needs a small amount of work; but it's not totally
        ## clear what uses it. I think that there's a good case for
        ## filling in nicely the requested bits - see sircovid's
        ## helper-lancelot-pmcmc.R which does this, and similar code
        ## in spimalot.
        if (predict$is_continuous) {
          times <- attr(predict$filter$data, "times")
          time <- c(times[[1]], times[, 2])
          trajectories <- mcstate_trajectories_continuous(
            time, trajectories_state, predicted = FALSE)
        } else {
          steps <- attr(predict$filter$data, "steps")
          step <- c(steps[[1]], steps[, 2])
          trajectories <- mcstate_trajectories_discrete(
            step, predict$rate, trajectories_state, predicted = FALSE)
        }
      }

      iteration <- seq(private$control$n_burnin + 1,
                       by = private$control$n_steps_every,
                       length.out = private$control$n_steps_retain)
      mcstate_pmcmc(iteration, pars, probabilities, state,
                    trajectories, restart, predict)
    }
  ))


history_collector <- function(n) {
  data <- vector("list", n)
  add <- function(i, value) {
    data[[i]] <<- value
  }

  get <- function() {
    data
  }

  list(add = add, get = get)
}


update_single <- function(f) {
  function(i) f()
}


update_alternate <- function(f, g, ratio) {
  if (ratio < 1) {
    return(update_alternate(g, f, 1 / ratio))
  }

  function(i) {
    if (i %% (ratio + 1) == 0) {
      g()
    } else {
      f()
    }
  }
}


make_rerun <- function(every, random) {
  if (!is.finite(every)) {
    function(i) FALSE
  } else if (random) {
    function(i) runif(1) < 1 / every
  } else {
    function(i) i %% every == 0
  }
}
