## In the case where we run different chains on different workers,
## this pair of classes takes care of the communication and scheduling
## details.
pmcmc_orchestrator <- R6::R6Class(
  "pmcmc_orchestrator",

  private = list(
    launch = function(process_id) {
      is_pending <- private$status == "pending"
      if (!any(is_pending)) {
        return()
      }
      chain_id <- which(is_pending)[[1]]
      private$target[[process_id]] <- chain_id
      private$status[[chain_id]] <- "running"
      private$sessions[[process_id]] <- callr::r_bg(
        pmcmc_chains_run,
        list(chain_id, private$path, private$n_threads[[chain_id]]),
        package = "mcstate")
    },

    control = NULL,
    path = NULL,

    sessions = NULL,  # list of session
    target = NULL,    # session -> chain map
    status = NULL,    # string, one of pending / running / done
    steps = NULL,     # number of steps complete
    n_threads = NULL, # number of threads per chain (precomputed)

    progress = NULL
  ),

  public = list(
    initialize = function(pars, initial, filter, control) {
      ## Changes behaviour of the internal progress bar so that we can
      ## interpret it as per-chain progress.
      control$progress_simple <- TRUE
      ## This is required, and we always set it
      control$use_parallel_seed <- TRUE
      ## We set this so that chains prepare doesn't get concerned;
      ## we've done the correct bookkeeping here and n_threads will be
      ## correct (but first compute the number of threads we'll use
      ## for chain)
      private$n_threads <- pmcmc_parallel_threads(
        control$n_threads_total, control$n_workers, control$n_chains)
      control$n_workers <- 1L

      path <- control$path %||% tempfile()
      private$path <- pmcmc_chains_prepare(path, pars, filter, control, initial)
      private$control <- control

      private$target <- rep(NA_integer_, control$n_chains)
      private$status <- rep("pending", control$n_chains)
      private$steps <- rep(0, control$n_chains)

      private$progress <-
        pmcmc_parallel_progress(control, private$status, private$steps)

      for (process_id in seq_len(control$n_workers)) {
        private$launch(process_id)
      }
    },

    ## processx::poll will poll, with a timeout, all our processes.
    ## There's not much downside to a long poll because if they *are*
    ## ready then they will return instantly. However, the process
    ## will only be interruptable each time the timeout triggers, so
    ## use 1000 here (1s).
    step = function(timeout = 1000) {
      res <- processx::poll(private$sessions, timeout)
      is_done <- vcapply(res, "[[", "process") == "ready"
      if (any(is_done)) {
        for (process_id in which(is_done)) {
          callr_safe_result(private$sessions[[process_id]])
          chain_id <- private$target[[process_id]]
          private$status[[chain_id]] <- "done"
          private$steps[[chain_id]] <- private$control$n_steps
          private$launch(process_id)
        }
      }

      if (private$control$progress) {
        has_stderr <- vcapply(res, "[[", "error") == "ready" & !is_done
        if (any(has_stderr)) {
          for (process_id in which(has_stderr)) {
            stderr <- private$sessions[[process_id]]$read_error_lines()
            re <- "^progress: ([0-9]+)$"
            i <- grep(re, stderr)
            if (length(i) > 0) {
              progress <- as.numeric(sub(re, "\\1", stderr[[last(i)]]))
              chain_id <- private$target[[process_id]]
              private$steps[[chain_id]] <- progress
            }
          }
        }
        private$progress(private$status, private$steps)
      }

      all(private$status == "done")
    },

    run = function() {
      while (!self$step()) {
      }
    },

    finish = function() {
      pmcmc_chains_collect(private$path)
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

  ## Grab another source of independent numbers to create the R
  ## seeds. This is essentially (though not identically) the behaviour
  ## of mcstate <= 0.6.16 which drew one number for the R seed from
  ## each generator but here we draw them all from the first.
  ##
  ## An alternative approach would be to take one long jump then
  ## generate seeds from that independent generator, but that has the
  ## downside of the R seed for each chain being dependent on the
  ## number of chains run.
  ##
  ## We rescale the real number to an integer on 1..2^24 (16777216)
  ## for the R seed as this will always integer representable.
  rng <- dust::dust_rng$new(seed_dust[[1]])$long_jump()
  seed_r <- ceiling(rng$random_real(n) * 2^24)

  Map(list, dust = seed_dust, r = seed_r)
}


## Create a callback to create a progress bar
pmcmc_parallel_progress <- function(control, status, steps, force = FALSE) {
  n_steps <- control$n_steps
  if (control$progress) {
    steps_total <- n_steps * control$n_chains
    fmt <- "[:spin] [:bar_overall] ETA :eta | :elapsedfull so far (:p_running)"
    t0 <- Sys.time()
    callback <- function(p) {
      message(sprintf("Finished %d steps in %s",
                      steps_total, format(Sys.time() - t0, digits = 1)))
    }
    p <- progress::progress_bar$new(fmt, steps_total, callback = callback,
                                    force = force)
    tick <- function(status, steps) {
      d <- pmcmc_parallel_progress_data(status, steps, n_steps)
      tryCatch(p$update(d$steps / steps_total, d$tokens),
               error = function(e) NULL)
    }

    ## Progress likes to be started right away:
    tick(status, steps)

    tick
  } else {
    function(status, steps) {
    }
  }
}

pmcmc_parallel_progress_data <- function(status, steps, n_steps) {
  map <- c(pending = " ", running = "+", done = "#")
  bar_overall <- paste(map[status], collapse = "")
  progress <- floor(steps[status == "running"] / n_steps * 100)
  p_running <- paste(sprintf("%3d%%", progress), collapse = " ")
  tokens <- list(bar_overall = bar_overall, p_running = p_running)
  list(steps = sum(steps), tokens = tokens)
}


## It is possible that we could allocate out threads better here
## really.  For the last wave of workers we would multiply n_threads
## by floor(n_threads * n_workers / n_workers_last); to get really
## fancy you could work out what the floor had discarded and add that
## to the last worker (which will be submitted last and benefit the
## most).
pmcmc_parallel_threads <- function(n_threads_total, n_workers, n_chains) {
  n_threads <- rep(n_threads_total / n_workers,
                   n_chains %/% n_workers * n_workers)
  spare <- n_chains %% n_workers
  if (spare != 0) {
    last <- c(rep(floor(n_threads_total / spare), spare - 1),
              floor(n_threads_total / spare) + n_threads_total %% spare)
    n_threads <- c(n_threads, last)
  }
  n_threads
}
