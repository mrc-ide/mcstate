## Things for working with our pmcmc objects

##' Combine multiple chains into a single (master) chain
##'
##' @title Combine multiple chains into a single (master) chain
##'
##' @param x An mcstate_pmcmc_list containing chains to combine
##'
##' @param burn_in an integer denoting the number of samples to discard
##' from each chain
##'
##' @export
create_master_chain <- function(x, burn_in) {

  if (class(x) != "mcstate_pmcmc_list") {
    stop("x must be a pmcmc_list object")
  }
  if (!is.numeric(burn_in)) {
    stop("burn_in must be an integer")
  }
  if (burn_in < 0) {
    stop("burn_in must not be negative")
  }
  if (burn_in >= nrow(x$chains[[1]]$results)) {
    stop("burn_in is greater than chain length")
  }

  chains <- lapply(
    X = x$chains,
    FUN = function(z) z$results[-seq_len(burn_in), ]
  )

  master <- do.call(what = rbind, args = chains)
  master
}

#
# Summary and plotting functions
#

##' @export
##' @importFrom stats cor sd
summary.mcstate_pmcmc <- function(object, ...) {

  par_names <- colnames(object$results)
  traces <- object$results

  # calculate correlation matrix
  corr_mat <- round(suppressWarnings(cor(traces)), 2)

  # compile summary
  summ <- rbind(mean = colMeans(traces),
                apply(traces, MARGIN = 2, quantile, c(0.025, 0.975)),
                min = apply(traces, MARGIN = 2, min),
                max =  apply(traces, MARGIN = 2, max)
  )
  summ <- as.data.frame(summ)
  summ <- round(summ, 3)

  sds <- round(apply(traces, 2, sd), 3)
  # convert start_date back into dates

  out <- list("summary" = summ,
              "corr_mat" = corr_mat,
              "sd" = sds)
  out
}

##' @export
summary.mcstate_pmcmc_list <- function(object, ..., burn_in = 1) {

  master_chain <- create_master_chain(x = object,
                                      burn_in = burn_in)

  z <- list(results = master_chain)
  summary.mcstate_pmcmc(z)
}

print_summ <- function(par_name, summ) {
  x <- summ$summary
  paste0(x["mean", par_name],
          "\n(",
          x["2.5%", par_name],
          ", ",
          x["97.5%", par_name], ")")
}

##' @export
##' @importFrom viridis cividis
##' @importFrom graphics hist par plot.new text
plot.mcstate_pmcmc <- function(x, ...) {

  summ <- summary(x)
  par_names <- colnames(x$results)
  n_pars <- length(par_names)

  traces <- x$results
  cols <- viridis::cividis(nrow(traces))
  cols <- cols[order(order(x$results$log_likelihood))]

  par(bty = "n",
      mfcol = c(n_pars, n_pars + 1L),
      mar = c(2.5, 2.5, 2, 1.5),
      mgp = c(1.5, 0.5, 0),
      oma = c(1, 1, 1, 1))

  for (i in seq_len(n_pars)) {
    for (j in seq_len(n_pars)) {
      if (i == j) {
        # plot hists on diagonal
        par_name <- par_names[i]
        breaks <- 10
        hist(traces[[i]],
             main = print_summ(par_name, summ),
             xlab = par_name,
             breaks = breaks,
             cex.main = 1,
             font.main = 1,
             freq = FALSE)
      } else if (i < j) {
        # plot correlations on lower triangle
        plot(x = traces[[i]],
             y = traces[[j]],
             xlab = par_names[i],
             ylab = par_names[j],
             col = cols,
             pch = 20)
      } else if (i > j) {
        # print rho on upper triangle
        plot.new()
        text(x = 0.5,
             y = 0.5,
             labels = paste("r =",
                            summ$corr_mat[i, j]))
      }
    }
  }

  # print traces in final column
  mapply(FUN = plot, traces,
         type = "l",
         ylab = par_names,
         xlab = "Iteration")
}

##' @export
##' @importFrom viridis cividis
##' @importFrom graphics hist par plot.new text lines legend matplot
plot.mcstate_pmcmc_list <- function(x, burn_in = 1, ...) {

  summ <- summary(x, burn_in = burn_in)
  par_names <- colnames(x$chains[[1]]$results)
  n_pars <- length(par_names)

  chains <- x$chains
  n_chains <- length(chains)
  cols_trace <- rev(viridis::viridis(n_chains))

  # compile master chain and order by log posterior for plotting
  master_chain <- create_master_chain(x, burn_in = burn_in)

  master_chain <- master_chain[order(master_chain$log_posterior), ]
  cols <- viridis::cividis(nrow(master_chain))
  cols <- cols[order(master_chain$log_posterior)]

  traces <- lapply(par_names, FUN = function(par_name) {
    lapply(X = chains,
           FUN = function(z) z$results[-seq_len(burn_in), par_name])
  })
  names(traces) <- par_names

  breaks <- lapply(par_names, function(par_name) {
    seq(from = min(master_chain[, par_name]),
        to =  max(master_chain[, par_name]),
        length.out = 20)
  })
  names(breaks) <- par_names

  hists <- lapply(par_names, FUN = function(par_name) {
    lapply(X = traces[[par_name]],
           FUN = hist,
           plot = FALSE,
           breaks = breaks[[par_name]])
  })
  names(hists) <- par_names

  hist_ylim <- lapply(hists, function(h) {
    chain_max <- sapply(h, function(chain) max(chain$density))
    upper_lim <- max(chain_max)
    if (is.na(upper_lim)) {
      upper_lim <- 0
    }
    c(0, upper_lim)
  })

  par(bty = "n",
      mfcol = c(n_pars, n_pars + 1L),
      mar = c(2.5, 2.5, 1.5, 0),
      mgp = c(1.5, 0.5, 0),
      oma = c(1, 1, 1, 1))

  for (i in seq_len(n_pars)) {
    for (j in seq_len(n_pars)) {
      if (i == j) {
        # plot hists on diagonal
        par_name <- par_names[i]
        bs <- breaks[[par_name]]
        plot(x = bs[1],
             y = 1,
             type = "n",
             xlim = c(bs[1], bs[length(bs)]),
             ylim = hist_ylim[[par_name]],
             xlab = par_name,
             ylab = "",
             main = print_summ(par_name, summ),
             cex.main = 1,
             font.main = 1
        )

        mapply(FUN = function(h, col) {
          plot_hists(h = h,
                     col = col,
                     breaks = bs)},
               h = hists[[par_name]],
               col = cols_trace)

      } else if (i < j) {
        # plot correlations on lower triangle
        plot(x = master_chain[[i]],
             y = master_chain[[j]],
             xlab = par_names[i],
             ylab = par_names[j],
             col = cols,
             pch = 20)
      } else if (i > j) {
        # print rho on upper triangle
        plot.new()
        text(x = 0.5,
             y = 0.5, cex = 1.5,
             labels = paste("r =",
                            summ$corr_mat[i, j]))
      }
    }
  }

  # print traces in final column
  n_iter <- nrow(master_chain) / n_chains

  mapply(FUN = function(par_name, leg) {

    trace <- do.call(cbind, traces[[par_name]])
    matplot(x = seq_len(nrow(trace)),
            y = trace,
            type = "l",
            col = cols_trace,
            lty = 1,
            xlab = "Iteration",
            ylab = par_name, )

    if (leg) {
      legend("top",
             ncol = n_chains,
             legend = paste("Chain", seq_len(n_chains)),
             fill = cols_trace,
             bty = "n")
    }
  },
  par_name = par_names,
  leg = c(TRUE, rep(FALSE, length(par_names) - 1)))
}

plot_hists <- function(h, col, breaks) {
  with(h, lines(x =  breaks,
                y = c(density,
                      density[length(density)]),
                type = "s",
                col = col))
}
