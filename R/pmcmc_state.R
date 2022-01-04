pmcmc_state <- R6::R6Class(
  "pmcmc_state",

  private = list(
    filter = NULL,
    pars = NULL,
    control = NULL,
    deterministic = NULL,

    nested = NULL,

    ## TODO: group together in a list?
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

    update_particle_history = function() {
      if (private$deterministic) {
        i <- 1L
      } else {
        i <- sample.int(private$filter$n_particles, 1)
      }

      ## TODO: there's some inconsistency throughout how we are
      ## accessing history; we should look into if it is expected that
      ## this will be dropped or not.

      if (private$control$save_trajectories) {
        private$curr_trajectories <- array_drop(private$filter$history(i), 2)
      }
      if (private$control$save_state) {
        private$curr_state <-
          array_drop(array_nth_dimension(private$filter$state(), 2, i), 2)
      }
      if (length(private$control$save_restart) > 0) {
        ## TODO: Can do this better by dropping within restart_state,
        ## which makes much more sense and is currently a bit weird.
        private$curr_restart <- array_drop(private$filter$restart_state(i), 2)
      }
    },

    update_mcmc_history = function() {
      private$history_pars$add(private$curr_pars)
      private$history_probabilities$add(
        c(private$curr_lprior, private$curr_llik, private$curr_lpost))

      if (!is.null(private$history_trajectories)) {
        private$history_trajectories$add(private$curr_trajectories)
      }
      if (!is.null(private$history_state)) {
        private$history_state$add(private$curr_state)
      }
      if (!is.null(private$history_restart)) {
        private$history_restart$add(private$curr_restart)
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
      private$curr_lpost - prop_lprior + log(u)
    },

    run_filter = function(p, min_log_likelihood = -Inf) {
      private$filter$run(private$pars$model(p),
                         private$control$save_trajectories,
                         private$control$save_restart,
                         min_log_likelihood)
    },

    update_simple = function() {
      prop_pars <- private$pars$propose(private$curr_pars)
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
    },

    update_fixed = function() {
      if (length(r6_private(private$pars)$fixed_parameters) > 0) {
        prop_pars <- private$pars$propose(private$curr_pars, type = "fixed")
        prop_lprior <- private$pars$prior(prop_pars)
        prop_llik <- private$run_filter(prop_pars)
        prop_lpost <- prop_lprior + prop_llik

        if (runif(1) < exp(sum(prop_lpost) - sum(private$curr_lpost))) {
          private$curr_pars <- prop_pars
          private$curr_lprior <- prop_lprior
          private$curr_llik <- prop_llik
          private$curr_lpost <- prop_lpost
          private$update_particle_history()
        }
      }
    },

    update_varied = function() {
      if (length(r6_private(private$pars)$varied_parameters) > 0) {
        prop_pars <- private$pars$propose(private$curr_pars, type = "varied")
        prop_lprior <- private$pars$prior(prop_pars)
        prop_llik <- private$run_filter(prop_pars)
        prop_lpost <- prop_lprior + prop_llik

        which <- runif(length(prop_lpost)) <
          exp(prop_lpost - private$curr_lpost)
        if (any(which)) {
          private$curr_pars[which, ] <- prop_pars[which, ]
          private$curr_lprior[which] <- prop_lprior[which]
          private$curr_llik[which] <- prop_llik[which]
          private$curr_lpost[which] <- prop_lpost[which]
          private$update_particle_history()
        }
      }
    },

    run_simple = function() {
      control <- private$control
      to <- min(private$curr_step + control$n_steps_each, control$n_steps)
      steps <- seq(from = private$curr_step + 1L,
                   length.out = to - private$curr_step)

      for (i in steps) {
        private$tick()

        if (rerun(i, control$rerun_every, control$rerun_random)) {
          private$curr_llik <- private$run_filter(private$curr_pars)
          private$curr_lpost <- private$curr_lprior + private$curr_llik
          private$update_particle_history()
        }

        private$update_simple()

        private$update_mcmc_history()
      }

      private$curr_step <- to

      list(step = to, finished = to == control$n_steps)
    },

    run_nested = function() {
      control <- private$control
      to <- min(private$curr_step + control$n_steps_each, control$n_steps)
      steps <- seq(from = private$curr_step + 1L,
                   length.out = to - private$curr_step)

      ## This bit is different, otherwise looking similar enough for
      ## dependency injection to work?
      run_alternate_step <- alternate(private$update_fixed,
                                      private$update_varied,
                                      control$nested_step_ratio)

      for (i in steps) {
        private$tick()

        if (rerun(i, control$rerun_every, control$rerun_random)) {
          private$curr_llik <- private$run_filter(private$curr_pars)
          private$curr_lpost <- private$curr_lprior + private$curr_llik
          private$update_particle_history()
        }

        run_alternate_step(i)

        private$update_mcmc_history()
      }
      private$curr_step <- to
      list(step = to, finished = to == control$n_steps)
    },

    finish_simple = function() {
      pars <- array_from_list(private$history_pars$get(), 2:1)
      colnames(pars) <- names(private$curr_pars)

      probabilities <- array_from_list(private$history_probabilities$get(), 2:1)
      colnames(probabilities) <-
        c("log_prior", "log_likelihood", "log_posterior")

      predict <- state <- restart <- trajectories <- NULL

      if (private$control$save_state || private$control$save_trajectories) {
        ## TODO: tidy up private access here; check what uses this?
        ##
        ## Do we *definitely* need step and rate here?
        data <- private$filter$inputs()$data
        predict <- list(transform = r6_private(private$pars)$transform,
                        index = r6_private(private$filter)$last_history$index,
                        step = last(data$step_end),
                        rate = attr(data, "rate", exact = TRUE),
                        filter = private$filter$inputs())
      }

      if (private$control$save_state) {
        state <- array_from_list(private$history_state$get(), 1:2)
      }

      if (length(private$control$save_restart) > 0) {
        ## [state x mcmc_sample x save point]
        restart_state <-
          array_from_list(private$history_restart$get(), c(1, 3, 2))
        restart <- list(time = private$control$save_restart,
                        state = restart_state)
      }

      if (private$control$save_trajectories) {
        ## [state x mcmc_sample x time]
        trajectories_state <-
          array_from_list(private$history_trajectories$get(), c(1, 3, 2))
        rownames(trajectories_state) <- names(predict$index)
        data <- private$filter$inputs()$data
        step <- c(data$step_start[[1]], data$step_end)
        trajectories <- mcstate_trajectories(step, predict$rate,
                                             trajectories_state, FALSE)
      }

      mcstate_pmcmc(pars, probabilities, state, trajectories, restart, predict)
    },

    ## TODO: this might be a good place to invert the structure of the
    ## generated parameters; we currently do so that each row is each
    ## pop, not each col
    finish_nested = function() {
      browser()
      ## pop x pop x step
      pars <- array_from_list(private$history_pars$get(), c(2, 1, 3))
      dimnames(pars)[2:1] <- dimnames(private$curr_pars)

      ## var x pop x step
      ## TODO: push this into the save
      tmp <- lapply(private$history_probabilities$get(), matrix,
                    ncol = 3, byrow = TRUE)
      probabilities <- array_from_list(tmp, c(2, 1, 3))
      dimnames(probabilities)[2:1] <-
        list(c("log_prior", "log_likelihood", "log_posterior"),
             rownames(private$curr_pars))

      predict <- state <- restart <- trajectories <- NULL

      if (private$control$save_state || private$control$save_trajectories) {
        data <- private$filter$inputs()$data
        predict <- list(transform = r6_private(private$pars)$transform,
                        index = r6_private(private$filter)$last_history$index,
                        step = last(data$step_end),
                        rate = attr(data, "rate", exact = TRUE),
                        filter = private$filter$inputs())
      }

      if (private$control$save_state) {
        # [state x pop x step]
        state <- list_to_array(private$history_state$get())
        colnames(state) <- rownames(private$curr_pars)
      }

      if (length(private$control$save_restart) > 0) {
        ## [state x sample x pop x time]
        restart_state  <-
          aperm(list_to_array(private$history_restart$get()), c(1, 4, 2, 3))
        restart <- list(time = private$control$save_restart,
                        state = restart_state)
      }

      if (private$control$save_trajectories) {
        ## [state x particle x population x step]
        trajectories_state <-
          list_to_array(private$history_trajectories$get())
        trajectories_state <- aperm(trajectories_state, c(1, 4, 2, 3))
        rownames(trajectories_state) <- names(predict$index)

        data <- private$filter$inputs()$data
        step_end <- data$step_end[seq_len(tabulate(data$population)[1])]
        step <- c(data$step_start[[1]], step_end)
        trajectories <- mcstate_trajectories(step, predict$rate,
                                             trajectories_state, FALSE)
      }

      mcstate_pmcmc(pars_array, probabilities, state, trajectories, restart,
                    predict)
    }
  ),

  public = list(
    initialize = function(pars, initial, filter, control) {
      private$filter <- filter
      private$pars <- pars
      private$control <- control
      private$nested <- inherits(pars, "pmcmc_parameters_nested")
      private$deterministic <- inherits(filter, "particle_deterministic")

      private$tick <- pmcmc_progress(control$n_steps, control$progress)

      private$curr_step <- 0L
      private$curr_pars <- initial
      private$curr_lprior <- private$pars$prior(private$curr_pars)
      private$curr_llik <- private$run_filter(private$curr_pars)
      private$curr_lpost <- private$curr_lprior + private$curr_llik
      private$update_particle_history()

      n_mcmc <- control$n_steps
      private$history_pars <- history_collector(n_mcmc)
      private$history_probabilities <- history_collector(n_mcmc)
      if (control$save_trajectories) {
        private$history_trajectories <- history_collector(n_mcmc)
      }
      if (control$save_state) {
        private$history_state <- history_collector(n_mcmc)
      }
      if (length(control$save_restart) > 0) {
        private$history_restart <- history_collector(n_mcmc)
      }

      private$update_mcmc_history()
    },

    set_n_threads = function(n_threads) {
      private$filter$set_n_threads(n_threads)
    },

    run = function() {
      if (private$nested) {
        private$run_nested()
      } else {
        private$run_simple()
      }
    },

    finish = function() {
      if (private$nested) {
        private$finish_nested()
      } else {
        private$finish_simple()
      }
    }
  ))


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

alternate <- function(f, g, ratio) {
  if (ratio < 1) {
    return(alternate(g, f, 1 / ratio))
  }

  function(i) {
    if (i %% (ratio + 1) == 0) {
      g()
    } else {
      f()
    }
  }
}


rerun <- function(i, every, random) {
  if (!is.finite(every)) {
    FALSE
  } else if (random) {
    runif(1) < 1 / every
  } else {
    i %% every == 0
  }
}
