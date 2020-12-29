orchestrator <- R6::R6Class(
  "orchestrator",

  private = list(
    n_steps = NULL,
    remote = NULL,
    sessions = NULL,
    status = NULL,
    results = NULL,
    n_threads = NULL,

    show_progress = function() {
      i <- lengths(private$status) > 0
      finished <- vlapply(private$status[i], "[[", "finished")
      n <- vnapply(private$status[i][!finished], "[[", "step")
      frac <- paste(sprintf("%d%%", round(n / private$n_steps * 100)),
                    collapse = " ")
      ## This is hard to read!
      message(sprintf("F %s R %s W %s | %s",
                      sum(finished), sum(!finished), sum(!i), frac))
      invisible(all(finished) && length(finished) == length(i))
    }
  ),

  public = list(
    ## TODO: some coordination will be required here to make sure that
    ## this aligns with the list in pmcmc() - at test at
    ## least.
    initialize = function(pars, filter, n_steps, save_state = TRUE,
                          save_trajectories = FALSE, progress = FALSE,
                          n_chains = 1, initial = NULL,
                          rerun_every = Inf,
                          ## Orchestrator arguments
                          n_workers = NULL, n_steps_each = NULL,
                          n_threads = NULL) {
      assert_is(pars, "pmcmc_parameters")
      assert_is(filter, "particle_filter")
      assert_scalar_positive_integer(n_steps)
      assert_scalar_positive_integer(n_steps_each)
      assert_scalar_logical(save_state)
      assert_scalar_logical(save_trajectories)
      assert_scalar_positive_integer(n_chains)

      if (n_chains < n_workers) {
        stop(sprintf("'n_chains' (%d) is less than 'n_workers' (%d)",
                     n_chains, n_workers))
      }

      ## Neither of these situations make any sense really, should we fail?
      ## if (n_chains < 2 || n_workers < 2) {
      ##   stop("What are you up to?")
      ## }

      inputs <- filter$inputs()
      seed <- make_seeds(n_chains, inputs$seed)
      initial <- pmcmc_check_initial(initial, pars, n_chains)

      if (is.null(n_threads)) {
        inputs$n_threads <- 1L
      } else {
        assert_scalar_positive_integer(n_threads)
        if (n_threads < n_workers) {
          stop(sprintf("'n_threads' (%d) is less than 'n_workers' (%d)",
                       n_threads, n_workers))
        }
        if (n_threads %% n_workers != 0) {
          stop(sprintf(
            "'n_threads' (%d) is not a multiple of 'n_workers' (%d)",
            n_threads, n_workers))
        }
        inputs$n_threads <- n_threads / n_workers
        private$n_threads <- list(n_threads = n_threads,
                                  n_workers = n_workers,
                                  free = 0L,
                                  target = n_threads / n_workers)
      }

      private$remote <- vector("list", n_workers)
      private$sessions <- vector("list", n_workers)
      private$status <- vector("list", n_chains)
      private$results <- vector("list", n_chains)
      control <- list(n_steps = n_steps,
                      n_steps_each = n_steps_each,
                      rerun_every = rerun_every,
                      save_state = save_state,
                      save_trajectories = save_trajectories)

      private$n_steps <- n_steps

      ## First stage starts the process, but this is async...
      for (i in seq_len(n_workers)) {
        private$remote[[i]] <- remote$new(pars, initial, inputs, control, seed)
        private$sessions[[i]] <- private$remote[[i]]$session
      }
      ## ...so once the sessions start coming up we start them working
      for (i in seq_len(n_workers)) {
        private$remote[[i]]$wait_session_ready()
        private$status[[i]] <- private$remote[[i]]$init(i)
      }
    },

    step = function() {
      res <- processx::poll(private$sessions, 1000)
      i <- vcapply(res, "[[", "process") == "ready"
      if (any(i)) {
        dat <- lapply(private$remote[i], function(x) x$read())
        index <- vnapply(dat, "[[", "index")
        result <- lapply(dat, function(x) x$data$result)
        private$status[index] <- result
        finished <- vlapply(result, function(x) x$finished)
        for (r in private$remote[i][!finished]) {
          private$n_threads <- update_threads_remove_from_pool(
            r, private$n_threads)
          r$continue()
        }
        if (any(finished)) {
          ## It might be preferable to close out sessions *first* as
          ## if we close out 2 sessions here simultaneously we
          ## underallocate cores for one time-chunk. That does
          ## complicate the book-keeping though.
          remaining <- which(lengths(private$status) == 0)
          for (r in private$remote[i][finished]) {
            res <- r$finish()
            private$results[res$index] <- list(res$data)
            if (length(remaining) == 0L) {
              r$session$close()
              private$n_threads <- update_threads_add_to_pool(
                r$n_threads, private$n_threads)
            } else {
              private$n_threads <- update_threads_remove_from_pool(
                r, private$n_threads)
              j <- remaining[[1]]
              remaining <- remaining[-1]
              private$status[[j]] <- r$init(j)
            }
          }
        }
      }
      private$show_progress()
    },

    get_results = function() {
      pmcmc_combine(samples = private$results)
    }
  ))


remote <- R6::R6Class(
  "remote",
  private = list(
    pars = NULL,
    initial = NULL,
    inputs = NULL,
    control = NULL,
    seed = NULL,
    step = NULL
  ),

  public = list(
    session = NULL,
    index = NULL,
    n_threads = NULL,

    initialize = function(pars, initial, inputs, control, seed) {
      self$session <- callr::r_session$new(wait = FALSE)

      private$pars <- pars
      private$initial <- initial
      private$inputs <- inputs
      private$control <- control
      private$seed <- seed

      lockBinding("session", self)
    },

    wait_session_ready = function(timeout = 3000) {
      self$session$poll_process(timeout)
      self$session$read()
      TRUE
    },

    ## Initialise this remote with the i'th chain (based on our
    ## initial conditions and seed that we began with). Every remote
    ## is *capable* of starting every chain but we do the allocation
    ## dynamically.
    init = function(index) {
      args <- list(private$pars, private$initial[, index], private$inputs,
                   private$control, private$seed[[index]])
      self$session$call(function(pars, initial, inputs, control, seed) {
        set.seed(seed$r)
        filter <- particle_filter_from_inputs(inputs, seed$dust)
        n_steps <- control$n_steps
        n_steps_each <- control$n_steps_each
        rerun_every <- control$rerun_every
        save_state <- control$save_state
        save_trajectories <- control$save_trajectories
        progress <- FALSE
        .GlobalEnv$obj <- pmcmc_state$new(
          pars, initial, filter, n_steps, n_steps_each, rerun_every,
          save_state, save_trajectories, progress)
        .GlobalEnv$obj$run()
      }, args, package = "mcstate")
      self$index <- index
      self$n_threads <- private$inputs$n_threads
      list(step = 0L, finished = FALSE)
    },

    continue = function() {
      self$session$call(function() .GlobalEnv$obj$run())
    },

    read = function() {
      list(index = self$index,
           data = self$session$read())
    },

    set_n_threads = function(n_threads) {
      message(sprintf("Setting to use %d threads", n_threads))
      self$session$run(function(n) .GlobalEnv$obj$set_n_threads(n),
                       list(n_threads))
      self$n_threads <- n_threads
    },

    ## This one is synchronous
    finish = function() {
      list(index = self$index,
           data = self$session$run(function() .GlobalEnv$obj$finish()))
    }
  ))


## We will, *before* starting anything, fully create a set of
## seeds. We use the dust rng to take a series of long_jumps (one per
## independent realisation), and from this also generate a single
## integer to use as the R seed.
##
## This means that the entire process is deterministic based on the
## single seed, which is itself sensibly chosen. Whether or not this
## leads to a sensible initialisation for R's RNG is a different
## question, but this should do as well as most reasonable attempts.
make_seeds <- function(n, seed) {
  ret <- vector("list", n)
  n_generators <- 1L
  if (is.raw(seed)) {
    n_generators <- length(seed) / 32L # 4 uint64_t, each 8 bytes
  }
  rng <- dust::dust_rng$new(seed, n_generators)
  for (i in seq_len(n)) {
    state <- rng$state()
    ret[[i]] <- list(dust = state,
                     r = rng$runif(1, 0, 1) * .Machine$integer.max)
    rng$long_jump()
  }
  ret
}


## This is the part when we drop workers off and return threads to the
## pool:
update_threads_add_to_pool <- function(n, dat) {
  if (is.null(dat)) {
    return(NULL)
  }
  dat$free <- dat$free + n
  dat$n_workers <- dat$n_workers - 1L
  dat$target <- ceiling(dat$n_threads / dat$n_workers)
  dat
}


update_threads_remove_from_pool <- function(r, dat) {
  if (is.null(dat) || dat$free == 0) {
    return(dat)
  }
  n <- min(r$n_threads + dat$free, dat$target)
  d <- n - r$n_threads
  if (d > 0) {
    r$set_n_threads(n)
    dat$free <- dat$free - d
  }

  dat
}
