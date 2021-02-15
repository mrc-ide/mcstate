pmcmc_state <- R6::R6Class(
  "pmcmc_state",

  private = list(
    filter = NULL,
    pars = NULL,
    control = NULL,

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

    update_history = function() {
      i <- sample.int(private$filter$n_particles, 1)
      if (private$control$save_trajectories) {
        private$curr_trajectories <-
          sample_trajectory(private$filter$history(i))
      }
      if (private$control$save_state) {
        private$curr_state <- private$filter$state()[, i, drop = TRUE]
      }
      if (length(private$control$save_restart) > 0) {
        private$curr_restart <- private$filter$restart_state(i)
        dim(private$curr_restart) <- dim(private$curr_restart)[-2L]
      }
    },

    run_filter = function(p) {
      private$filter$run(private$pars$model(p),
                         private$control$save_trajectories,
                         private$control$save_restart)
    }
  ),

  public = list(
    initialize = function(pars, initial, filter, control) {
      private$filter <- filter
      private$pars <- pars
      private$control <- control

      private$history_pars <- history_collector(control$n_steps)
      private$history_probabilities <- history_collector(control$n_steps)
      private$history_state <- history_collector(control$n_steps)
      private$history_restart <- history_collector(control$n_steps)
      private$history_trajectories <- history_collector(control$n_steps)

      private$curr_step <- 0L
      private$curr_pars <- initial
      private$curr_lprior <- private$pars$prior(private$curr_pars)
      private$curr_llik <- private$run_filter(private$curr_pars)
      private$curr_lpost <- private$curr_lprior + private$curr_llik

      private$history_pars$add(private$curr_pars)
      private$history_probabilities$add(c(private$curr_lprior,
                                          private$curr_llik,
                                          private$curr_lpost))

      private$tick <- pmcmc_progress(control$n_steps, control$progress)

      ## Initial version of the history
      private$update_history()
      if (private$control$save_trajectories) {
        private$history_trajectories$add(private$curr_trajectories)
      }
      if (private$control$save_state) {
        private$history_state$add(private$curr_state)
      }
      if (length(private$control$save_restart) > 0) {
        private$history_restart$add(private$curr_restart)
      }
    },

    set_n_threads = function(n_threads) {
      private$filter$set_n_threads(n_threads)
    },

    run = function() {
      to <- min(private$curr_step + private$control$n_steps_each,
                private$control$n_steps)
      steps <- seq(from = private$curr_step + 1L,
                   length.out = to - private$curr_step)
      for (i in steps) {
        private$tick()

        if (i %% private$control$rerun_every == 0) {
          private$curr_llik <- private$run_filter(private$curr_pars)
          private$curr_lpost <- private$curr_lprior + private$curr_llik
          private$update_history()
        }

        prop_pars <- private$pars$propose(private$curr_pars)
        prop_lprior <- private$pars$prior(prop_pars)
        prop_llik <- private$run_filter(prop_pars)
        prop_lpost <- prop_lprior + prop_llik

        if (runif(1) < exp(prop_lpost - private$curr_lpost)) {
          private$curr_pars <- prop_pars
          private$curr_lprior <- prop_lprior
          private$curr_llik <- prop_llik
          private$curr_lpost <- prop_lpost
          private$update_history()
        }

        private$history_pars$add(private$curr_pars)
        private$history_probabilities$add(
          c(private$curr_lprior, private$curr_llik, private$curr_lpost))
        if (private$control$save_trajectories) {
          private$history_trajectories$add(private$curr_trajectories)
        }
        if (private$control$save_state) {
          private$history_state$add(private$curr_state)
        }
        if (length(private$control$save_restart) > 0) {
          private$history_restart$add(private$curr_restart)
        }
      }
      private$curr_step <- to
      list(step = to, finished = to == private$control$n_steps)
    },

    finish = function() {
      pars_matrix <- set_colnames(list_to_matrix(private$history_pars$get()),
                                  names(private$curr_pars))
      probabilities <- set_colnames(
        list_to_matrix(private$history_probabilities$get()),
        c("log_prior", "log_likelihood", "log_posterior"))

      predict <- state <- restart <- trajectories <- NULL

      if (private$control$save_state || private$control$save_trajectories) {
        ## Do we *definitely* need step and rate here?
        data <- private$filter$inputs()$data
        predict <- list(transform = r6_private(private$pars)$transform,
                        index = r6_private(private$filter)$last_history$index,
                        step = last(data$step_end),
                        rate = attr(data, "rate", exact = TRUE),
                        filter = private$filter$inputs())
      }

      if (private$control$save_state) {
        state <- t(list_to_matrix(private$history_state$get()))
      }

      if (length(private$control$save_restart) > 0) {
        ## Permute state from [state x save point x mcmc_sample]
        ## to [state x mcmc_sample x save point] to match the predicted state)
        restart_state  <-
          aperm(list_to_array(private$history_restart$get()), c(1, 3, 2))
        restart <- list(time = private$control$save_restart,
                        state = restart_state)
      }

      if (private$control$save_trajectories) {
        ## Permute trajectories from [state x mcmc x particle] to
        ## [state x particle x mcmc] so that they match the ones that we
        ## will generate with predict
        trajectories_state <-
          aperm(list_to_array(private$history_trajectories$get()), c(1, 3, 2))
        rownames(trajectories_state) <- names(predict$index)
        data <- private$filter$inputs()$data
        step <- c(data$step_start[[1]], data$step_end)
        trajectories <- mcstate_trajectories(step, predict$rate,
                                             trajectories_state, FALSE)
      }

      mcstate_pmcmc(pars_matrix, probabilities, state, trajectories, restart,
                    predict)
    }
  ))


## A utility function for sampling a trajectory and safely dropping
## the dimensionality even if there is only one state vector
sample_trajectory <- function(history, index) {
  ret <- history[, index, , drop = TRUE]
  if (is.null(dim(ret))) {
    dim(ret) <- dim(history)[c(1, 3)]
  }
  ret
}


## Generic history collector, collects anything at all into a list
##
## This would be more nicely done as a simple R6 class but it's a bit
## slow in testing; this version speeds up the total mcmc runtime by a
## factor of ~3x (0.4s/1000 iterations to 0.13s/1000) mostly by
## reducing the number of garbage collections considerably.
history_collector <- function(n) {
  data <- vector("list", n + 1L)
  i <- 0L
  add <- function(value) {
    i <<- i + 1L
    data[[i]] <<- value
  }

  get <- function() {
    data
  }

  list(add = add, get = get)
}