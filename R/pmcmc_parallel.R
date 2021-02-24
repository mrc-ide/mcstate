## In the case where we run different chains on different workers,
## this pair of classes takes care of the communication and scheduling
## details.
pmcmc_orchestrator <- R6::R6Class(
  "pmcmc_orchestrator",

  private = list(
    control = NULL,
    remotes = NULL,
    sessions = NULL,
    status = NULL,
    results = NULL,
    thread_pool = NULL,
    progress = NULL,
    nested = FALSE
  ),

  public = list(
    initialize = function(pars, initial, filter, control) {
      private$control <- control
      private$thread_pool <- thread_pool$new(control$n_threads_total,
                                             control$n_workers)
      private$progress <- pmcmc_parallel_progress(control)

      inputs <- filter$inputs()
      seed <- make_seeds(control$n_chains, inputs$seed)
      inputs$n_threads <- private$thread_pool$target

      private$remotes <- vector("list", control$n_workers)
      private$sessions <- vector("list", control$n_workers)
      private$status <- vector("list", control$n_chains)
      private$results <- vector("list", control$n_chains)
      if (inherits(pars, "pmcmc_parameters_nested")) {
        private$nested <- TRUE
      }

      ## First stage starts the process, but this is async...
      for (i in seq_len(control$n_workers)) {
        private$remotes[[i]] <- pmcmc_remote$new(
          pars, initial, inputs, control, seed)
        private$sessions[[i]] <- private$remotes[[i]]$session
      }
      ## ...so once the sessions start coming up we start them working
      for (i in seq_len(control$n_workers)) {
        private$remotes[[i]]$wait_session_ready()
        private$status[[i]] <- private$remotes[[i]]$init(i)
      }
    },

    step = function() {
      ## processx::poll will poll, with a timeout, our
      ## processes. There's not much downside to a long poll because
      ## if they *are* ready they will return instantly. However, the
      ## process will only be interruptable each time the timeout
      ## triggers, so use 1000 here (1s).
      res <- processx::poll(private$sessions, 1000)
      i <- vcapply(res, "[[", "process") == "ready"
      if (any(i)) {
        dat <- lapply(private$remotes[i], function(x) x$read())
        index <- vnapply(dat, "[[", "index")
        result <- lapply(dat, "[[", "result")
        private$status[index] <- result
        finished <- vlapply(result, function(x) x$finished)
        for (r in private$remotes[i][!finished]) {
          private$thread_pool$remove(r)
          r$continue()
        }
        if (any(finished)) {
          ## It might be preferable to close out sessions *first* as
          ## if we close out 2 sessions here simultaneously we
          ## underallocate cores for one time-chunk. That does
          ## complicate the book-keeping though.
          remaining <- which(lengths(private$status) == 0)
          for (r in private$remotes[i][finished]) {
            res <- r$finish()
            private$results[res$index] <- list(res$data)
            if (length(remaining) == 0L) {
              r$session$close()
              private$thread_pool$add(r)
            } else {
              private$thread_pool$remove(r)
              j <- remaining[[1]]
              remaining <- remaining[-1]
              private$status[[j]] <- r$init(j)
            }
          }
        }
      }
      private$progress(private$status)
    },

    run = function() {
      while (!self$step()) {
      }
    },

    finish = function() {
      if (private$nested) {
        pmcmc_combine_nested(samples = private$results)
      } else {
        pmcmc_combine(samples = private$results)
      }
    }
  ))


## This class takes care of the details if a partially run chain
## running in a remote process, wrapping around callr's "r_session"
## objects.
pmcmc_remote <- R6::R6Class(
  "pmcmc_remote",
  private = list(
    pars = NULL,
    initial = NULL,
    inputs = NULL,
    control = NULL,
    seed = NULL,
    step = NULL,
    nested = FALSE
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
      if (inherits(pars, "pmcmc_parameters_nested")) {
        private$nested <- TRUE
      }

      lockBinding("session", self)
    },

    ## 3000ms is the timeout for the session to come alive; an error
    ## will be thrown if it is not alive by then (this number is used
    ## by callr internally)
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
      if (is_3d_array(private$initial)) {
        initial <- private$initial[, , index]
      } else {
        initial <- private$initial[, index]
      }
      args <- list(private$pars, initial, private$inputs,
                   private$control, private$seed[[index]])
      self$session$call(function(pars, initial, inputs, control, seed) {
        set.seed(seed$r)
        filter <- particle_filter_from_inputs(inputs, seed$dust)
        control$progress <- FALSE
        .GlobalEnv$obj <- pmcmc_state$new(pars, initial, filter, control)
        if (inherits(pars, "pmcmc_parameters_nested")) {
          .GlobalEnv$obj$run_nested()
        } else {
          .GlobalEnv$obj$run()
        }
      }, args, package = "mcstate")
      self$index <- index
      self$n_threads <- private$inputs$n_threads
      list(step = 0L, finished = FALSE)
    },

    continue = function() {
      if (private$nested) {
        self$session$call(function() .GlobalEnv$obj$run_nested())
      } else {
        self$session$call(function() .GlobalEnv$obj$run())
      }
    },

    read = function() {
      data <- self$session$read()
      if (!is.null(data$error)) {
        ## NOTE: We have to use this non-exported function to get the
        ## same nice error handling as Gabor has set up in the
        ## package, and depending on any of the details of the
        ## returned objects will likely be even more fragile than
        ## grabbing this function from within the package.
        callr:::throw(data$error)
      }
      list(index = self$index, result = data$result)
    },

    set_n_threads = function(n_threads) {
      self$session$run(function(n) .GlobalEnv$obj$set_n_threads(n),
                       list(n_threads))
      self$n_threads <- n_threads
    },

    ## This one is synchronous
    finish = function() {
      if (private$nested) {
        list(index = self$index,
             data = self$session$run(function()
               .GlobalEnv$obj$finish_nested()))
      } else {
        list(index = self$index,
             data = self$session$run(function()
               .GlobalEnv$obj$finish()))
      }
    }
  ))


## We will, *before* starting anything, fully create a set of seeds,
## one per chain regardless of how many workers are being used. We use
## the dust rng to take a series of long_jumps (one per independent
## realisation), and from this also generate a single integer to use
## as the R seed.
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


## Utility class to help with the thread book-keeping above.
thread_pool <- R6::R6Class(
  class = FALSE,
  public = list(
    n_threads = NULL,
    n_workers = NULL,
    free = 0L,
    target = NULL,
    active = NULL,

    initialize = function(n_threads, n_workers) {
      self$active <- !is.null(n_threads)
      self$n_threads <- n_threads
      self$n_workers <- n_workers
      self$free <- 0L

      if (self$active) {
        self$target <- ceiling(n_threads / n_workers)
      } else {
        self$target <- 1L
      }
    },

    add = function(remote) {
      if (!self$active) {
        return()
      }
      self$n_workers <- self$n_workers - 1L
      self$free <- self$free + remote$n_threads
      self$target <- ceiling(self$n_threads / self$n_workers)
    },

    remove = function(remote) {
      if (!self$active || self$free == 0L) {
        return()
      }
      n <- min(remote$n_threads + self$free, self$target)
      d <- n - remote$n_threads
      if (d > 0) {
        remote$set_n_threads(n)
        self$free <- self$free - d
      }
    }
  ))


pmcmc_parallel_progress_data <- function(status, n_steps) {
  started <- lengths(status) > 0
  finished <- vlapply(status[started], "[[", "finished")

  n_finished <- sum(finished)
  n_started <- sum(started)
  n_running <- n_started - n_finished
  n_waiting <- length(status) - n_started

  n <- vnapply(status[started], "[[", "step")

  ## Could use cli::symbol$full_block and crayon here to make this nicer
  bar_overall <- paste0(strrep("#", n_finished),
                        strrep("+", n_running),
                        strrep(" ", n_waiting))
  p_running <- paste(sprintf("%d%%", round(n[!finished] / n_steps * 100)),
                     collapse = " ")

  tokens <- list(bar_overall = bar_overall, p_running = p_running)
  result <- n_finished == length(status)
  list(n = sum(n), tokens = tokens, result = result)
}


## Create a callback to create a progress bar
pmcmc_parallel_progress <- function(control, force = FALSE) {
  n_steps <- control$n_steps
  if (control$progress) {
    n <- n_steps * control$n_chains
    fmt <- "[:spin] [:bar_overall] ETA :eta | :elapsedfull so far (:p_running)"
    t0 <- Sys.time()
    callback <- function(p) {
      message(sprintf("Finished %d steps in %s",
                      n, format(Sys.time() - t0, digits = 0)))
    }
    p <- progress::progress_bar$new(fmt, n, callback = callback, force = force)
    ## Progress likes to be started right away:
    d <- pmcmc_parallel_progress_data(vector("list", control$n_chains), n_steps)
    p$update(d$n, d$tokens)
    function(status) {
      d <- pmcmc_parallel_progress_data(status, n_steps)
      p$update(d$n / n, d$tokens)
      d$result
    }
  } else {
    function(status) {
      pmcmc_parallel_progress_data(status, n_steps)$result
    }
  }
}
