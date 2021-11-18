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
    path = NULL,
    filter_inputs = NULL,
    nested = FALSE
  ),

  public = list(
    initialize = function(pars, initial, filter, control, root = NULL) {
      private$control <- control
      private$thread_pool <- thread_pool$new(control$n_threads_total,
                                             control$n_workers)
      private$progress <- pmcmc_parallel_progress(control)

      filter_inputs <- filter$inputs()
      seed <- make_seeds(control$n_chains, filter_inputs$seed, filter$model)
      filter_inputs$n_threads <- private$thread_pool$target

      private$remotes <- vector("list", control$n_workers)
      private$sessions <- vector("list", control$n_workers)
      private$status <- vector("list", control$n_chains)
      private$results <- vector("list", control$n_chains)
      private$filter_inputs <- filter_inputs
      private$nested <- inherits(pars, "pmcmc_parameters_nested")

      root <- root %||% tempfile()
      dir.create(root, FALSE, TRUE)
      private$path <- list(root = root,
                           input = file.path(root, "input.rds"),
                           output = file.path(root, "output-%d.rds"))

      input <- list(
        pars = pars,
        initial = pmcmc_parallel_initial(control$n_chains, initial),
        filter = filter_inputs,
        control = control,
        seed = seed)

      ## Ignore warning:
      ##   'package:mcstate' may not be available when loading
      ## which would cause significantly more issues than here :)
      suppressWarnings(saveRDS(input, private$path$input))

      ## First stage starts the process and reads in input data, but
      ## this is async over the workers
      n_threads <- filter_inputs$n_threads
      nested <- inherits(pars, "pmcmc_parameters_nested")
      for (i in seq_len(control$n_workers)) {
        private$remotes[[i]] <-
          pmcmc_remote$new(private$path$input, n_threads, nested)
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
            filename <- sprintf(private$path$output, r$index)
            res <- r$finish(filename)
            dat <- readRDS(filename)
            private$results[[r$index]] <-
              pmcmc_parallel_predict_filter(dat, private$filter_inputs)
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
      unlink(private$path$root, recursive = TRUE)
      pmcmc_combine(samples = private$results)
    }
  ))


## This class takes care of the details if a partially run chain
## running in a remote process, wrapping around callr's "r_session"
## objects.
pmcmc_remote <- R6::R6Class(
  "pmcmc_remote",
  private = list(
    path = NULL,
    step = NULL,
    nested = NULL
  ),

  public = list(
    session = NULL,
    index = NULL,
    n_threads = NULL,

    ## NOTE: n_threads here must match that of the filter inputs
    initialize = function(path, n_threads, nested) {
      options <- callr::r_session_options(
        load_hook = bquote(.GlobalEnv$input <- readRDS(.(path))))
      self$session <- callr::r_session$new(options = options, wait = FALSE)
      self$n_threads <- n_threads
      private$nested <- nested
      lockBinding("session", self)
    },

    ## 30000ms is the timeout for the session to come alive; an error
    ## will be thrown if it is not alive by then (this number is used
    ## by callr internally)
    wait_session_ready = function(timeout = 30000) {
      self$session$poll_process(timeout)
      self$session$read()
      TRUE
    },

    ## Initialise this remote with the i'th chain (based on our
    ## initial conditions and seed that we began with). Every remote
    ## is *capable* of starting every chain but we do the allocation
    ## dynamically.
    init = function(index) {
      self$session$call(function(index, nested) {
        ## simplify resolution, technically not needed
        input <- .GlobalEnv$input
        seed <- input$seed[[index]]
        initial <- input$initial[[index]]
        control <- input$control

        set.seed(seed$r)
        filter <- particle_filter_from_inputs(input$filter, seed$dust)
        control$progress <- FALSE
        .GlobalEnv$obj <- pmcmc_state$new(input$pars, initial, filter, control)
        if (nested) {
          .GlobalEnv$obj$run_nested()
        } else {
          .GlobalEnv$obj$run()
        }
      }, list(index, private$nested), package = "mcstate")
      self$index <- index
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

    ## This one is synchronous, and writes to disk. Using callr's I/O
    ## here is too slow. We might want to make this async, but it will
    ## really complicate the above!
    finish = function(filename) {
      method <- if (private$nested) "finish_nested" else "finish"

      self$session$run(function(method, filename) {
        results <- .GlobalEnv$obj[[method]]()
        results$predict$filter <- results$predict$filter$seed
        suppressWarnings(saveRDS(results, filename))
      }, list(method, filename))

      list(index = self$index, data = filename)
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
make_seeds <- function(n, seed, model) {
  n_streams <- 1L
  if (is.raw(seed)) {
    ## This is not always correct and varies with the model; move to
    ## make this explicit I think in the next step.
    n_streams <- length(seed) / 32L # 4 uint64_t, each 8 bytes
  }

  seed_dust <- dust::dust_rng_distributed_state(seed, n_streams, n, model)

  ## Grab another source of independent numbers to create the R seeds
  ## with by doing one further long jump from the last state
  ##
  ## Using integers on 1..2^24 (16777216) for the R seed as always
  ## integer representable.
  rng <- dust::dust_rng$new(seed_dust[[1]])$long_jump()
  seed_r <- ceiling(rng$random_real(n) * 2^24)

  Map(list, dust = seed_dust, r = seed_r)
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
                      n, format(Sys.time() - t0, digits = 1)))
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


pmcmc_parallel_initial <- function(n_chains, initial) {
  if (is_3d_array(initial)) {
    initial <- lapply(seq_len(n_chains), function(index)
      initial[, , index])
  } else {
    initial <- lapply(seq_len(n_chains), function(index)
      initial[, index])
  }
  initial
}


pmcmc_parallel_predict_filter <- function(dat, filter_inputs) {
  if (!is.null(dat$predict$filter)) {
    filter_inputs$seed <- dat$predict$filter
    dat$predict$filter <- filter_inputs
  }
  dat
}
