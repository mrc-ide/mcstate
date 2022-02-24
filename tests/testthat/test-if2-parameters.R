context("IF2 (parameters)")

test_that("can construct a parameter", {
  p <- if2_parameter("a", 0.1, min = 0, max = 1)

  expect_s3_class(p, "if2_parameter")
  expect_equal(p$name, "a")
  expect_equal(p$min, 0)
  expect_equal(p$max, 1)
  expect_equal(p$initial, 0.1)
  expect_false(p$integer)
  expect_equal(p$prior(), 0)
})


test_that("if2_parameter must satisfy min/max constraints", {
  expect_error(
    p <- if2_parameter("a", -1, min = 0, max = 1),
    "'initial' must be >= 'min' (0)",
    fixed = TRUE)
  expect_error(
    p <- if2_parameter("a", 2, min = 0, max = 1),
    "'initial' must be <= 'max' (1)",
    fixed = TRUE)
})


test_that("if2_parameter prior works", {
  expect_error(
    p <- if2_parameter("a", 0, min = -1, max = 1,
                       prior = function(x) 1 / x),
    "Prior function for 'a' returned a non-finite value on initial value",
    fixed = TRUE)
  expect_error(
    p <- if2_parameter("a", -1, min = -1, max = 1,
                       prior = function(x) sample(c(0, 1), x)),
    paste0("Prior function for 'a' failed to evaluate initial value:",
           " invalid 'size' argument"),
    fixed = TRUE)
})


test_that("can construct and walk a set of parameters", {
  pars <- if2_parameters$new(
            list(if2_parameter("beta", 0.15, min = 0.1, max = 2),
                 if2_parameter("gamma", 0.05, min = 0, max = 1),
                 if2_parameter("time", 10, integer = TRUE)))
  expect_s3_class(pars, "if2_parameters")
  expect_equal(pars$names(), c("beta", "gamma", "time"))
  expect_equal(pars$initial(), c("beta" = 0.15, "gamma" = 0.05, "time" = 10))
  expect_equal(pars$summary(),
               data_frame(name = c("beta", "gamma", "time"),
                          min = c(0.1, 0, -Inf),
                          max = c(2, 1, Inf),
                          integer = c(FALSE, FALSE, TRUE)))

  n_pars <- length(pars$names())
  n_par_sets <- 5
  pars_sd <- c(0.1, 0.1, 2)
  walk_mat <- pars$walk_initialise(n_par_sets, pars_sd)
  expect_equal(dim(walk_mat), c(n_pars, n_par_sets))
  expect_true(all(walk_mat[1, ] >= 0.1 & walk_mat[1, ] <= 2))
  expect_true(all(walk_mat[2, ] >= 0 & walk_mat[2, ] <= 1))
  # Check rounded correctly
  expect_true(all(vlapply(walk_mat[3, ],
              function(x) {
                max(abs(round(x) - x)) < sqrt(.Machine$double.eps)
              })))

  model_input <- pars$model(walk_mat)
  expect_equal(length(model_input), n_par_sets)
  expect_equal(length(model_input[[1]]), n_pars)
})
