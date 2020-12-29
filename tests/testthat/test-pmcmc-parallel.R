context("pmmc (parallel)")

test_that("basic parallel operation", {
  dat <- example_sir()
  n_particles <- 42
  n_steps <- 30

  p0 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = dat$index, seed = 1L)
  s <- make_seeds(2L, 1L)

  obj <- orchestrator$new(dat$pars, p0, n_steps, n_chains = 2,
                          n_workers = 2, n_steps_each = 20)
  while (!obj$step()) {
  }
  ans <- obj$get_results()

  ## Run two chains manually with a given pair of seeds:
  s <- make_seeds(2L, 1L)
  set.seed(s[[1]]$r)
  p1 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = dat$index, seed = s[[1]]$dust)
  res1 <- pmcmc(dat$pars, p1, n_steps)
  set.seed(s[[2]]$r)
  p2 <- particle_filter$new(dat$data, dat$model, n_particles, dat$compare,
                            index = dat$index, seed = s[[2]]$dust)
  res2 <- pmcmc(dat$pars, p2, n_steps)

  ## Combined set:
  res3 <- pmcmc_combine(res1, res2)

  expect_equal(res3$pars, ans$pars)
  expect_equal(res3, ans)
})


test_that("make seeds with integer returns length-32 raws", {
  s <- make_seeds(2, 1)
  expect_length(s, 2L)
  expect_setequal(names(s[[1]]), c("dust", "r"))
  expect_length(s[[1]]$dust, 32L)
  expect_is(s[[1]]$dust, "raw")
  expect_identical(make_seeds(2, 1), s)
  expect_identical(make_seeds(4, 1)[1:2], s)
})


test_that("make seeds with long raw retains size", {
  seed <- dust::dust_rng$new(NULL, n_generators = 3)$state()
  s <- make_seeds(2L, seed)
  expect_length(s, 2L)
  expect_setequal(names(s[[1]]), c("dust", "r"))
  expect_identical(s[[1]]$dust, seed)
  expect_length(s[[2]]$dust, length(seed))
})
