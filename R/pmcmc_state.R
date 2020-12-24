pmcmc_state <- R6::R6Class(
  "pmcmc_state",

  private = list(
    filter = NULL,
    pars = NULL,
    rerun_every = NULL,
    n_steps = NULL,
    n_particles = NULL,
    save_state = NULL,
    save_trajectories = NULL,

    history_pars = NULL,
    history_probabilities = NULL,
    history_state = NULL,
    history_trajectories = NULL,

    curr_pars = NULL,
    curr_lprior = NULL,
    curr_llik = NULL,
    curr_lpost = NULL,
    curr_trajectories = NULL,
    curr_state = NULL,

    tick = NULL,

    update_history = function() {
      i <- sample.int(private$n_particles, 1)
      if (private$save_trajectories) {
        private$curr_trajectories <-
          sample_trajectory(private$filter$history(i))
      }
      if (private$save_state) {
        private$curr_state <- private$filter$state()[, i, drop = TRUE]
      }
    },

    run_filter = function(p) {
      private$filter$run(private$pars$model(p), private$save_trajectories)
    }
  ),

  public = list(
    initialize = function(pars, initial, filter, n_steps, rerun_every,
                          save_state, save_trajectories, progress) {
      private$filter <- filter
      private$pars <- pars

      private$n_steps <- n_steps
      private$n_particles <- filter$n_particles
      private$rerun_every <- rerun_every
      private$save_state <- save_state
      private$save_trajectories <- save_trajectories

      private$history_pars <- history_collector(n_steps)
      private$history_probabilities <- history_collector(n_steps)
      private$history_state <- history_collector(n_steps)
      private$history_trajectories <- history_collector(n_steps)

      private$curr_pars <- initial
      private$curr_lprior <- private$pars$prior(private$curr_pars)
      private$curr_llik <- filter$run(private$pars$model(private$curr_pars),
                                      private$save_trajectories)
      private$curr_lpost <- private$curr_lprior + private$curr_llik

      private$history_pars$add(private$curr_pars)
      private$history_probabilities$add(c(private$curr_lprior,
                                          private$curr_llik,
                                          private$curr_lpost))

      private$tick <- pmcmc_progress(n_steps, progress)

      ## Initial version of the history
      private$update_history()
      if (private$save_trajectories) {
        private$history_trajectories$add(private$curr_trajectories)
      }
      if (private$save_state) {
        private$history_state$add(private$curr_state)
      }
    },

    ## This needs to become easier to run partially
    run = function() {
      for (i in seq_len(private$n_steps)) {
        private$tick()

        if (i %% private$rerun_every == 0) {
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
        if (private$save_trajectories) {
          private$history_trajectories$add(private$curr_trajectories)
        }
        if (private$save_state) {
          private$history_state$add(private$curr_state)
        }
      }

      pars_matrix <- set_colnames(list_to_matrix(private$history_pars$get()),
                                  names(private$curr_pars))
      probabilities <- set_colnames(
        list_to_matrix(private$history_probabilities$get()),
        c("log_prior", "log_likelihood", "log_posterior"))

      predict <- state <- trajectories <- NULL

      if (private$save_state || private$save_trajectories) {
        ## Do we *definitely* need step and rate here?
        data <- private$filter$inputs()$data
        predict <- list(transform = r6_private(private$pars)$transform,
                        index = r6_private(private$filter)$last_index_state,
                        step = last(data$step_end),
                        rate = attr(data, "rate", exact = TRUE),
                        filter = private$filter$inputs())
      }

      if (private$save_state) {
        state <- t(list_to_matrix(private$history_state$get()))
      }

      if (private$save_trajectories) {
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

      mcstate_pmcmc(pars_matrix, probabilities, state, trajectories, predict)
    }
  ))
