## Each node in the tree has elements:
##
## * parent - the index of node from which this node comes from
## * theta - the index of the node to take parameters from if we fail
## * p - the probability of this node being reached, given p_accept
## * children - index of child nodes on reject and accept (respectively)
##
## The tree will be constructed such that the probabilities 'p' are
## decreasing (with ties)
speculate_tree <- function(n, p_accept) {
  tree <- vector("list", n)

  node <- function(parent, theta, p) {
    list(parent = parent,
         theta = theta,
         p = p,
         children = c(NA_integer_, NA_integer_))
  }

  children <- function(parent, tree) {
    list(
      list(p = tree[[parent]]$p * (1 - p_accept),
           accept = FALSE,
           parent = parent),
      list(p = tree[[parent]]$p * p_accept,
           accept = TRUE,
           parent = parent))
  }

  tree[[1]] <- node(0, 0, 1)
  live <- children(1, tree)
  live_p <- c(live[[1]]$p, live[[2]]$p)

  for (i in seq_len(n)[-1]) {
    j <- which.max(live_p)
    x <- live[[j]]
    tree[[x$parent]]$children[[x$accept + 1L]] <- i
    theta <- if (x$accept) x$parent else tree[[x$parent]]$theta
    tree[[i]] <- node(x$parent, theta, x$p)
    kids <- children(i, tree)
    live <- c(live[-j], kids)
    live_p <- c(live_p[-j], c(kids[[1]]$p, kids[[2]]$p))
  }

  tree
}


speculate_propose <- function(theta, proposal, tree, include_current) {
  n <- length(tree)
  ret <- matrix(NA_real_, n, length(theta))
  for (i in seq_len(n)) {
    x <- tree[[i]]
    ret[i, ] <- proposal(if (x$theta == 0) theta else ret[x$parent, ])
  }
  ret
}


speculate_accept <- function(theta, p, theta_new, p_new, tree) {
  chain_index <- integer()
  chain_accept <- logical()
  curr <- 0
  idx <- 1L
  while (!is.na(idx)) {
    accept <- p_new[[idx]] > p || runif(1) < exp(p_new[[idx]] - p)
    if (accept) {
      p <- p_new[[idx]]
      curr <- idx
    }
    chain_index <- c(chain_index, curr)
    chain_accept <- c(chain_accept, accept)
    idx <- tree[[idx]]$children[[accept + 1]]
  }

  list(length = length(chain_index),
       index = chain_index,
       accept = chain_accept)
}


speculate_tree_truncate <- function(tree) {
  n <- length(tree)
  ret <- tree[-n]
  parent <- ret[[tree[[n]]$parent]]
  parent$children[parent$children == n] <- NA_integer_
  ret[[tree[[n]]$parent]] <- parent
  ret
}


speculate_prepare_filter <- function(filter, control) {
  if (filter$has_multiple_parameters) {
    if (filter$n_parameters != control$speculate_n) {
      stop("Your filter is configured for an incompatible number of parameters")
    }
  } else {
    inputs <- filter$inputs()
    inputs$n_parameters <- control$speculate_n
    filter <- particle_filter_from_inputs(inputs)
  }
  filter
}
