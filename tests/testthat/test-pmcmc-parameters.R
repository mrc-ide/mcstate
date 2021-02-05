context("pmcmc_parameters")

test_that("Can construct a parameter", {
  p <- pmcmc_parameter("a", 1, 0, 10)
  expect_s3_class(p, "pmcmc_parameter")
  expect_equal(p$name, "a")
  expect_equal(p$initial, 1)
  expect_equal(p$min, 0)
  expect_equal(p$max, 10)
  expect_false(p$discrete)
  expect_equal(p$prior(1), 0)
})


test_that("Can provide a prior", {
  f <- function(p) log(1 / p)
  p <- pmcmc_parameter("a", 1, prior = f)
  expect_identical(p$prior, f)
})


test_that("parameter initial point must lie in range", {
  expect_silent(pmcmc_parameter("a", 0))
  expect_silent(pmcmc_parameter("a", 0, min = 0))
  expect_silent(pmcmc_parameter("a", 0, max = 0))
  expect_error(pmcmc_parameter("a", 0, min = 1),
               "'initial' must be >= 'min' (1)", fixed = TRUE)
  expect_error(pmcmc_parameter("a", 0, max = -1),
               "'initial' must be <= 'max' (-1)", fixed = TRUE)
})


test_that("initial value must satify prior and not fail", {
  expect_silent(
    pmcmc_parameter("a", 10, prior = log))
  expect_error(
    pmcmc_parameter("a", 0, prior = log),
    "Prior function for 'a' returned a non-finite value on initial value")
  expect_error(
    pmcmc_parameter("a", 0, prior = function(p) stop("an error")),
    "Prior function for 'a' failed to evaluate initial value: an error")
})


test_that("can scale proposal kernels", {
  proposal_kernel <- diag(2) * 1e-4
  row.names(proposal_kernel) <- colnames(proposal_kernel) <- c("beta", "gamma")

  p <- pmcmc_parameters$new(
    list(
      pmcmc_parameter("beta", 0.2, min = 0, max = 1,
                      prior = function(p) log(1e-10)),
      pmcmc_parameter("gamma", 0.1, min = 0, max = 1,
                      prior = function(p) log(1e-10))),
    proposal = proposal_kernel)

  expect_equal(
    p$initial(),
    c(beta = 0.2, gamma = 0.1))

  set.seed(1)
  p0 <- p$initial()
  p1 <- p$propose(p0)
  set.seed(1)
  p2 <- p$propose(p0, 10)

  set.seed(1)
  expect_equal(p1, rmvnorm_generator(proposal_kernel)(p0))

  expect_equal(p2 - p0, (p1 - p0) * sqrt(10))
})


test_that("transform parameters returns a list by default", {
  proposal_kernel <- diag(2) * 1e-4
  row.names(proposal_kernel) <- colnames(proposal_kernel) <- c("beta", "gamma")

  p <- pmcmc_parameters$new(
    list(
      pmcmc_parameter("beta", 0.2, min = 0, max = 1,
                      prior = function(p) log(1e-10)),
      pmcmc_parameter("gamma", 0.1, min = 0, max = 1,
                      prior = function(p) log(1e-10))),
    proposal = proposal_kernel)

  expect_equal(
    p$model(p$initial()),
    list(beta = 0.2, gamma = 0.1))
  expect_equal(
    p$model(p$initial() + 0.5),
    list(beta = 0.7, gamma = 0.6))
})


test_that("provide custom transform and use it", {
  proposal_kernel <- diag(2) * 1e-4
  row.names(proposal_kernel) <- colnames(proposal_kernel) <- c("beta", "gamma")

  transform <- function(p) {
    list(alpha = p[["beta"]] + p[["gamma"]],
         gamma = rep(p[["gamma"]], 4),
         beta = p[["beta"]])
  }

  p <- pmcmc_parameters$new(
    list(pmcmc_parameter("beta", 0.2), pmcmc_parameter("gamma", 0.1)),
    proposal = proposal_kernel,
    transform = transform)

  expect_equal(
    p$model(p$initial()),
    list(alpha = 0.3, gamma = rep(0.1, 4), beta = 0.2))
  expect_equal(
    p$model(p$initial() + 0.5),
    list(alpha = 1.3, gamma = rep(0.6, 4), beta = 0.7))
})


test_that("can compute prior", {
  proposal_kernel <- diag(2) * 1e-4
  row.names(proposal_kernel) <- colnames(proposal_kernel) <- c("beta", "gamma")

  transform <- function(p) {
    list(alpha = p[["beta"]] + p[["gamma"]],
         gamma = rep(p[["gamma"]], 4),
         beta = p[["beta"]])
  }

  p <- pmcmc_parameters$new(
    list(
      pmcmc_parameter("beta", 0.2, prior = function(p) dexp(p, log = TRUE)),
      pmcmc_parameter("gamma", 0.1, prior = function(p) dnorm(p, log = TRUE))),
    proposal = diag(2),
    transform = transform)

  ## Compute prior
  expect_equal(
    p$prior(p$initial()),
    dexp(0.2, log = TRUE) + dnorm(0.1, log = TRUE))
  ## Done by position, not name:
  expect_equal(
    p$prior(list(0.2, 0.1)),
    dexp(0.2, log = TRUE) + dnorm(0.1, log = TRUE))
  ## Non-default values
  expect_equal(
    p$prior(list(beta = 1, gamma = 0)),
    dexp(1, log = TRUE) + dnorm(0, log = TRUE))
})


test_that("all inputs to pmcmc_parameters must be pmcmc_parameter objects", {
  expect_error(
    pmcmc_parameters$new(NULL, proposal = diag(2)),
    "'parameters' must be a list")
  expect_error(
    pmcmc_parameters$new(list(), proposal = diag(2)),
    "At least one parameter is required")

  expect_error(
    pmcmc_parameters$new(list(a = 1, b = 2), proposal = diag(2)),
    "Expected all elements of '...' to be 'pmcmc_parameter' objects",
    fixed = TRUE)
})


test_that("proposal matrix and parameters must conform", {
  pars <- list(pmcmc_parameter("beta", 0.2), pmcmc_parameter("gamma", 0.1))
  expect_error(
    pmcmc_parameters$new(pars, proposal = diag(1)),
    "Expected a square proposal matrix with 2 rows and columns")
  expect_error(
    pmcmc_parameters$new(pars, proposal = diag(3)),
    "Expected a square proposal matrix with 2 rows and columns")
  expect_error(
    pmcmc_parameters$new(pars, proposal = matrix(0, 3, 2)),
    "Expected a square proposal matrix with 2 rows and columns")
})


test_that("if proposal has names, they must match parameters", {
  pars <- list(pmcmc_parameter("beta", 0.2), pmcmc_parameter("gamma", 0.1))
  nms <- c("beta", "gamma")
  m <- function(a, b = a) {
    matrix(1, length(a), length(b), dimnames = list(a, b))
  }

  expect_error(
    pmcmc_parameters$new(pars, proposal = m(c("a", "b"))),
    "Expected dimension names of 'proposal' to match parameters")
  expect_error(
    pmcmc_parameters$new(pars, proposal = m(c("a", "b"), nms)),
    "Expected dimension names of 'proposal' to match parameters")
  expect_error(
    pmcmc_parameters$new(pars, proposal = m(nms, c("a", "b"))),
    "Expected dimension names of 'proposal' to match parameters")
  expect_silent(pmcmc_parameters$new(pars, proposal = m(nms)))
})


test_that("can summarise parameters", {
  p <- pmcmc_parameters$new(
    list(
      pmcmc_parameter("beta", 0.2, min = 0, max = 1,
                      prior = function(p) log(1e-10)),
      pmcmc_parameter("gamma", 0.1, min = 0, max = 1,
                      prior = function(p) log(1e-10))),
    proposal = diag(2))

  expect_equal(p$names(), c("beta", "gamma"))
  expect_equal(
    p$summary(),
    data_frame(name = c("beta", "gamma"), min = 0, max = 1, discrete = FALSE))
})


test_that("duplicate parameter names are rejected", {
  expect_error(
    pmcmc_parameters$new(
      list(pmcmc_parameter("a", 0.2), pmcmc_parameter("a", 0.1)),
      matrix(1, 2, 2)),
    "Duplicate parameter names: 'a'")
  expect_error(
    pmcmc_parameters$new(
      list(pmcmc_parameter("a", 0.2),
           pmcmc_parameter("a", 0.1),
           pmcmc_parameter("a", 0.1),
           pmcmc_parameter("b", 0.1),
           pmcmc_parameter("c", 0.1),
           pmcmc_parameter("c", 0.1)),
      matrix(1, 6, 6)),
    "Duplicate parameter names: 'a', 'c'")
})


test_that("named parameters must match parameter names", {
  pars <- list(pmcmc_parameter("a", 0.2), pmcmc_parameter("b", 0.1))
  m <- matrix(1, 2, 2)
  expect_silent(
    pmcmc_parameters$new(setNames(pars, c("a", "b")), m))
  expect_error(
    pmcmc_parameters$new(setNames(pars, c("b", "a")), m),
    "'parameters' is named, but the names do not match parameters")
})


test_that("can fix parameters", {
  n <- 5
  nms <- letters[seq_len(n)]
  vcv <- var(matrix(runif(n * n), n, n))
  rownames(vcv) <- colnames(vcv) <- nms

  initial <- runif(n)
  prior_rate <- runif(n)
  pars <- Map(function(nm, i, r)
    pmcmc_parameter(nm, i, prior = function(p) log(r)),
    nms, initial, prior_rate)

  p <- pmcmc_parameters$new(pars, proposal = vcv)
  p2 <- p$fix(c(b = 0.5, d = 0.2))

  expect_equal(p2$names(), c("a", "c", "e"))
  expect_equal(p2$initial(), p$initial()[c("a", "c", "e")])
  expect_equal(p2$propose(p2$initial(), 0), p2$initial())

  ## Loose check here; the underlying implementation is simple enough
  ## though
  i <- c(1, 3, 5)
  expect_equal(
    var(t(replicate(5000, p2$propose(initial[i])))),
    unname(vcv[i, i]),
    tolerance = 1e-2)

  cmp <- as.list(p$initial())
  cmp[c("b", "d")] <- c(0.5, 0.2)
  expect_equal(p2$model(p2$initial()), cmp)
})


test_that("prevent impossible fixed parameters", {
  n <- 5
  nms <- letters[seq_len(n)]
  vcv <- var(matrix(runif(n * n), n, n))
  rownames(vcv) <- colnames(vcv) <- nms

  initial <- runif(n)
  prior_rate <- runif(n)
  pars <- Map(function(nm, i, r)
    pmcmc_parameter(nm, i, prior = function(p) log(r)),
    nms, initial, prior_rate)

  p <- pmcmc_parameters$new(pars, proposal = vcv)
  expect_error(p$fix(c(1, 2)), "'fixed' must be named")
  expect_error(p$fix(c(a = 1, b = 2, a = 1)), "'fixed' must have unique names")
  expect_error(p$fix(c(a = 1, b = 2, f = 1)),
               "Fixed parameters not found in model: 'f'")
  expect_error(p$fix(c(a = 1, b = 1, c = 1, d = 1, e = 1)),
               "Cannot fix all parameters")
})

test_that("Can construct a varied parameter", {
  p <- pmcmc_varied_parameter("a", 1:10, 0, 10)
  expect_s3_class(p, "pmcmc_varied_parameter")
  expect_equal(p$name, "a")
  expect_equal(p$initial, 1:10)
  expect_equal(p$min, numeric(10))
  expect_equal(p$max, rep(10, 10))
  expect_false(p$discrete)
  expect_equal(p$prior[[1]](1), 0)
})

test_that("can construct varied parameters", {
 expect_silent(pmcmc_varied_parameter(name = "a", 1))
 expect_silent(pmcmc_varied_parameter(name = "a", 1:100))

  expect_silent(pmcmc_varied_parameter(name = "a", initial = 0.1, min = 0,
    max = 1, prior = function(a) log(a)))

  expect_silent(pmcmc_varied_parameter(name = "a", initial = 1:5, min = 0,
    max = 6, prior = function(a) log(a)))

  expect_error(pmcmc_varied_parameter(name = "a", initial = 1:5, min = 0,
    max = 4, prior = function(a) log(a)), "<= 'max'")

  expect_error(pmcmc_varied_parameter(name = "a", initial = 1:5, min = 2,
    max = 6, prior = function(a) log(a)), ">= 'min'")

  expect_error(pmcmc_varied_parameter(name = "a", initial = 1:5, min = 0:1,
    max = 6, prior = function(a) log(a)), "Length of 'min'")

  expect_error(pmcmc_varied_parameter(name = "a", initial = 1:5, min = 0,
    max = 3:4, prior = function(a) log(a)), "Length of 'max'")

  expect_equal(
    pmcmc_varied_parameter(name = "a", initial = 1:5, min = 0,
      max = 6, prior = list(function(a) log(a))),
    pmcmc_varied_parameter(name = "a", initial = 1:5, min = 0,
      max = 6, prior = function(a) log(a)))

  expect_error(pmcmc_varied_parameter(name = "a", initial = 1:5, min = 0,
    max = 6, prior = function(a) Inf), "non-finite value")

  expect_error(pmcmc_varied_parameter(name = "a", initial = 1:5, min = 0,
    max = 6, prior = function(a) stop()), "Prior function for 'a'")
})

test_that("varied parameters in pmcmc_parameters", {
  expect_silent(pmcmc_parameters$new(list(
    pmcmc_varied_parameter(name = "a", 1),
    pmcmc_parameter(name = "b", 1)
  ), matrix(c(1, 0.5, 0.5, 1), 2, 2), populations = c("Europe", "America")))

  expect_error(
    pmcmc_parameters$new(list(
      pmcmc_varied_parameter(name = "a", 1:3, 2:4, min = c(0, 1, 2)),
      pmcmc_parameter(name = "b", 1)
    ), matrix(c(1, 0.5, 0.5, 1), 2, 2), populations = c("Europe", "America")),
    "Expected length"
  )
})

test_that("varied pmcmc_parameters summary", {
  p <- pmcmc_parameters$new(
    parameters = list(
      pmcmc_varied_parameter(name = "a", 2:4, min = 0:2, max = 4:6),
      pmcmc_parameter(name = "b", 2, min = 1, max = 3)),
    proposal = matrix(c(1, 0.5, 0.5, 1), 2, 2),
    populations = c("Europe", "America", "Africa")
  )

  expect_identical(p$summary("Europe"),
    data.frame(name = c("a", "b"), min = c(0, 1), max = c(4, 3),
               discrete = FALSE, type = c("varied", "fixed"))
  )
  expect_identical(p$summary("Africa"),
    data.frame(name = c("a", "b"), min = c(2, 1), max = c(6, 3),
               discrete = FALSE, type = c("varied", "fixed"))
  )
  expect_error(p$summary("Australia"), "Expected 'population'")

  s <- p$summary()
  expect_equal(s,
              list(Europe = p$summary("Europe"),
                  America = p$summary("America"),
                  Africa = p$summary("Africa")))

})

test_that("varied initial", {
    p <- pmcmc_parameters$new(
    parameters = list(
      pmcmc_varied_parameter(name = "a", 2:4, min = 0:2, max = 4:6),
      pmcmc_parameter(name = "b", 2, min = 1, max = 3)),
    proposal = matrix(c(1, 0.5, 0.5, 1), 2, 2),
    populations = c("Europe", "America", "Africa")
  )
  expect_identical(p$initial(),
    matrix(c(2, 2, 3, 2, 4, 2), 3, 2, TRUE,
    dimnames = list(c("Europe", "America", "Africa"), letters[1:2])))
})


test_that("varied proposal", {
  p <- pmcmc_parameters$new(
    parameters = list(
      pmcmc_varied_parameter(name = "a", 2:4, min = 0:2, max = 4:6),
      pmcmc_parameter(name = "b", 2, min = 1, max = 3)),
    proposal = matrix(c(1, 0.5, 0.5, 1), 2, 2),
    populations = c("Europe", "America", "Africa")
  )
  init <- p$initial()

  aprop <- p$propose(init, type = "fixed")
  expect_equal(dimnames(aprop),
    list(c("Europe", "America", "Africa"), c("a", "b")))
  expect_equal(aprop[, 1], init[, 1])
  expect_true(all(aprop[, 2] != init[, 2]))

  aprop <- p$propose(init, type = "varied")
  expect_equal(dimnames(aprop),
    list(c("Europe", "America", "Africa"), c("a", "b")))
  expect_equal(aprop[, 2], init[, 2])
  expect_true(all(aprop[, 1] != init[, 1]))

  aprop <- p$propose(init, type = "both")
  expect_equal(dimnames(aprop),
    list(c("Europe", "America", "Africa"), c("a", "b")))
  expect_true(all(aprop[, 1] != init[, 1]))
  expect_true(all(aprop[, 2] != init[, 2]))
})

test_that("varied prior", {
    p <- pmcmc_parameters$new(
    parameters = list(
      pmcmc_varied_parameter(name = "a", 2:4, min = 0:2, max = 4:6,
        prior = list(function(x) 1, function(x) 2, function(x) 3)),
      pmcmc_parameter(name = "b", 2, min = 1, max = 3, prior = function(x) 4)),
    proposal = matrix(c(1, 0.5, 0.5, 1), 2, 2),
    populations = c("Europe", "America", "Africa")
  )
  init <- p$initial()
  expect_equal(p$prior(init), set_names(5:7, c("Europe", "America", "Africa")))

  p <- pmcmc_parameters$new(
    parameters = list(
      pmcmc_varied_parameter(name = "a", 2:4, min = 0:2, max = 4:6,
        prior = list(dnorm, dexp, function(x) dnorm(x, 2))),
      pmcmc_parameter(name = "b", 2, min = 1, max = 3, prior = dlnorm)),
    proposal = matrix(c(1, 0.5, 0.5, 1), 2, 2),
    populations = c("Europe", "America", "Africa")
  )
  init <- p$initial()
  expect_equal(p$prior(init),
    set_names(c(dnorm(2) + dlnorm(2), dexp(3) + dlnorm(2), dnorm(4, 2) + dlnorm(2)),
    c("Europe", "America", "Africa")))
})

test_that("varied model/transform", {
  p <- pmcmc_parameters$new(
    parameters = list(
      pmcmc_varied_parameter(name = "a", 2:4, min = 0:2, max = 4:6,
        prior = list(dnorm, dexp, function(x) dnorm(x, 2))),
      pmcmc_parameter(name = "b", 2, min = 1, max = 3, prior = dlnorm)),
    proposal = matrix(c(1, 0.5, 0.5, 1), 2, 2),
    populations = c("Europe", "America", "Africa")
  )
  expect_equal(p$model(p$initial()), apply(p$initial(), 1, as.list))

  p <- pmcmc_parameters$new(
    parameters = list(
      pmcmc_varied_parameter(name = "a", 2:4, min = 0:2, max = 4:6,
        prior = list(dnorm, dexp, function(x) dnorm(x, 2))),
      pmcmc_parameter(name = "b", 2, min = 1, max = 3, prior = dlnorm)),
    proposal = matrix(c(1, 0.5, 0.5, 1), 2, 2),
    populations = c("Europe", "America", "Africa"),
    transform = function(x) as.list(log(x))
  )
  expect_equal(p$model(p$initial()),
                apply(p$initial(), 1, function(x) as.list(log(x))))
})

test_that("varied fixed", {
  p <- pmcmc_parameters$new(
    parameters = list(
      pmcmc_varied_parameter(name = "a", 2:4, min = 0:2, max = 4:6,
        prior = list(dnorm, dexp, function(x) dnorm(x, 2))),
      pmcmc_parameter(name = "b", 2, min = 1, max = 3, prior = dlnorm),
      pmcmc_parameter(name = "c", 3),
      pmcmc_varied_parameter(name = "d", 4)),
    proposal = diag(4),
    populations = c("Europe", "America", "Africa")
  )
  p_fix <- pmcmc_parameters$new(
    parameters = list(
      pmcmc_parameter(name = "b", 2, min = 1, max = 3, prior = dlnorm),
      pmcmc_varied_parameter(name = "d", 4)),
    proposal = diag(2),
    populations = c("Europe", "America", "Africa"),
    transform = function(p) {
          base_transform(set_into(base, idx_vary, p))
      }
  )
  fixedp <- p$fix(c(a = 1, c = 2))
  expect_equal(fixedp, p_fix)

  expect_equal(
    fixedp$initial(),
    matrix(c(2, 4), nrow = 3, ncol = 2, TRUE,
    list(c("Europe", "America", "Africa"), c("b", "d")))
  )

  mod <- fixedp$model(fixedp$initial())
  expect_true(is.list(mod) && length(mod) == 3)
  expect_equal(names(mod), c("Europe", "America", "Africa"))
  expect_equal(mod$Europe, list(a = 1, b = 2, c = 2, d = 4))
})
