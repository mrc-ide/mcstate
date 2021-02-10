
test_that("Can construct a varied parameter", {
  p <- pmcmc_varied_parameter("p1", letters[1:2], 1:2)
  expect_s3_class(p, "pmcmc_varied_parameter")
  expect_equal(p$a$name, "p1")
  expect_equal(p$b$name, "p1")
  expect_equal(p$a$initial, 1)
  expect_equal(p$b$initial, 2)
  expect_equal(p$a$min, -Inf)
  expect_equal(p$b$min, - Inf)
  expect_equal(p$a$max, Inf)
  expect_equal(p$b$max, Inf)
  expect_equal(p$a$discrete, FALSE)
  expect_equal(p$b$discrete, FALSE)
  expect_equal(p$a$prior, function(p) 0)
  expect_equal(p$b$prior, function(p) 0)
})

test_that("varied parameter reps", {
  expect_equal(pmcmc_varied_parameter("p1", letters[1:3], 1)$c$initial, 1)
  expect_equal(
    pmcmc_varied_parameter("p1", letters[1:3], 1, min = 1)$c$min, 1)
  expect_equal(
    pmcmc_varied_parameter("p1", letters[1:3], 1, max = 1)$c$max, 1)
  expect_true(
    pmcmc_varied_parameter("p1", letters[1:3], 1, discrete = TRUE)$c$discrete)
  expect_equal(
    pmcmc_varied_parameter("p1", letters[1:3], 1,
        prior = function(x) 1)$c$prior,
    function(x) 1)
  expect_equal(
    pmcmc_varied_parameter("p1", letters[1:3], 1,
        prior = list(function(x) 1))$c$prior,
    function(x) 1)
})

test_that("varied parameter errors", {
  expect_error(pmcmc_varied_parameter("p1", letters[1:3], 1:2), "Length of")
  expect_error(pmcmc_varied_parameter("p1", letters[1:3], 1, min = 1:2),
    "Length of")
  expect_error(pmcmc_varied_parameter("p1", letters[1:3], 1, max = 1:2),
    "Length of")
  expect_error(
    pmcmc_varied_parameter("p1", letters[1:3], 1, discrete = c(TRUE, FALSE)),
    "Length of")
  expect_error(
    pmcmc_varied_parameter("p1", letters[1:3], 1,
        prior = rep(list(function(x) 1), 2)),
    "Length of")
})

test_that("clean_parameters", {
   expect_error(clean_parameters(list()), "At least one")
   expect_error(clean_parameters(list("a")), "Expected all")
   expect_error(clean_parameters(list(pmcmc_parameter("a", 2))),
       "No varied")
   expect_error(clean_parameters(
       list(pmcmc_varied_parameter("a", "a", 2),
            pmcmc_varied_parameter("a", "a", 3))),
         "Duplicate")
   expect_error(clean_parameters(
       list(pmcmc_varied_parameter("a", "a", 2),
            pmcmc_parameter("a", 3))),
         "Duplicate")
    expect_error(clean_parameters(
       list(a = pmcmc_varied_parameter("a", "a", 2),
            c = pmcmc_parameter("b", 3))),
         "Fixed parameters are named")
    expect_error(clean_parameters(
       list(c = pmcmc_varied_parameter("a", "a", 2),
            b = pmcmc_parameter("b", 3))),
         "Varied parameters are named")
    expect_error(clean_parameters(
       list(pmcmc_varied_parameter("a", c("p1", "p2"), 2),
            pmcmc_varied_parameter("b", c("p2", "p1"), 2),
            pmcmc_parameter("c", 3))),
         "Populations and ordering")
    expect_error(clean_parameters(
       list(pmcmc_varied_parameter("a", c("p1", "p2"), 2),
            pmcmc_varied_parameter("b", c("p2"), 2),
            pmcmc_parameter("c", 3))),
         "Populations and ordering")

    expect_silent(clean_parameters(
       list(a = pmcmc_varied_parameter("a", c("p1", "p2"), 2),
            b = pmcmc_varied_parameter("b", c("p1", "p2"), 2),
            c = pmcmc_parameter("c", 3))))
    expect_equal(clean_parameters(
       list(a = pmcmc_varied_parameter("a", c("p1", "p2"), 2),
        b = pmcmc_varied_parameter("b", c("p1", "p2"), 2)))$fixed_parameters,
        NULL)
})

test_that("clean proposals", {
    params <- list(fixed = letters[1:2], varied = letters[3:4])
    pops <- letters[5:6]

    expect_equal(
      clean_proposals(diag(2), diag(2), pops, params)$
        varied,
      array(diag(2), dim = c(2, 2, 2),
      dimnames = list(letters[3:4], letters[3:4], letters[5:6])))

     # edge case, 1 of each
    expect_equal(
      clean_proposals(diag(1), diag(1), pops, list(fixed = "a", varied = "b"))$
        varied,
      array(diag(1), dim = c(1, 1, 2),
      dimnames = list("b", "b", pops)))

    expect_equal(
      clean_proposals(diag(2), diag(2), pops, params)$
        fixed,
      matrix(diag(2), 2, 2, dimnames = list(letters[1:2], letters[1:2])))

    expect_error(
      clean_proposals(diag(2),
        array(diag(2), c(1, 2, 2),
            dimnames = list(NULL, NULL, letters[3:4])),
            pops, params),
      "Expected proposal")

    expect_error(
      clean_proposals(diag(2),
        array(diag(2), c(2, 2, 2),
            dimnames = list(NULL, NULL, letters[3:4])),
            pops, params),
      "Expected 3rd")

    expect_error(
      clean_proposals(diag(2),
        array(diag(2), c(2, 2, 2),
            dimnames = list(letters[1:2], letters[3:4], letters[5:6])),
            pops, params),
      "proposal_varied")

    expect_error(
      clean_proposals(
        matrix(diag(2), 2, 2, dimnames = list(letters[5:6], letters[5:6])),
        array(diag(2), c(2, 2, 2)), pops, params),
      "proposal_fixed")
})

test_that("construct pmcmc_nested_parameters", {
    parameters <- list(a = pmcmc_varied_parameter("a", c("p1", "p2"), 2),
                       b = pmcmc_varied_parameter("b", c("p1", "p2"), 2),
                       c = pmcmc_parameter("c", 3),
                       d = pmcmc_parameter("d", 4))
    proposal_fixed <- diag(2)
    proposal_varied <- diag(2) + 1

    expect_silent(pmcmc_nested_parameters$new(parameters, proposal_varied,
      proposal_fixed))
    expect_silent(pmcmc_nested_parameters$new(parameters, proposal_varied,
      proposal_fixed, transform = as.list))

    # edge case, 1 each
    parameters <- list(a = pmcmc_varied_parameter("a", c("p1", "p2"), 2),
                       d = pmcmc_parameter("d", 4))
    proposal_fixed <- diag(1)
    proposal_varied <- diag(1) + 1
    expect_silent(pmcmc_nested_parameters$new(parameters, proposal_varied,
      proposal_fixed))

    # varied only
    parameters <- list(a = pmcmc_varied_parameter("a", c("p1", "p2"), 2))
    proposal_varied <- diag(1) + 1
    expect_silent(pmcmc_nested_parameters$new(parameters, proposal_varied))
})

test_that("pmcmc_nested_parameters initial", {
  parameters <- list(a = pmcmc_varied_parameter("a", c("p1", "p2"), 1:2),
                     b = pmcmc_varied_parameter("b", c("p1", "p2"), 3:4),
                     c = pmcmc_parameter("c", 5),
                     d = pmcmc_parameter("d", 6))
  proposal_fixed <- diag(2)
  proposal_varied <- diag(2) + 1
  p <- pmcmc_nested_parameters$new(parameters, proposal_varied, proposal_fixed)
  expect_equal(
    p$initial(),
    matrix(c(5, 5, 6, 6, 1:4), nrow = 2, ncol = 4,
        dimnames = list(c("p1", "p2"), c("c", "d", "a", "b")))
  )

  # varied only
  parameters <- list(
    a = pmcmc_varied_parameter("a", c("p1", "p2"), 1:2),
    b = pmcmc_varied_parameter("b", c("p1", "p2"), 3:4)
  )
  proposal_varied <- diag(2)
  p <- pmcmc_nested_parameters$new(parameters, proposal_varied)
  expect_equal(
    p$initial(),
    matrix(1:4, 2, 2, dimnames = list(c("p1", "p2"), letters[1:2]))
  )

  # 1 varied 1 pop
  parameters <- list(
    a = pmcmc_varied_parameter("a", c("p1"), 1)
  )
  proposal_varied <- diag(1)
  p <- pmcmc_nested_parameters$new(parameters, proposal_varied)
  expect_equal(
    p$initial(),
    matrix(1, nrow = 1, ncol = 1, dimnames = list(c("p1"), "a"))
  )

  # 1 varied 1 fix
  parameters <- list(
    a = pmcmc_varied_parameter("a", c("p1", "p2"), 1, prior = dnorm),
    b = pmcmc_parameter("b", 1, prior = dexp))
  proposal_varied <- diag(1)
  proposal_fixed <- diag(1)
  p <- pmcmc_nested_parameters$new(parameters, proposal_varied, proposal_fixed)
  expect_equal(
    p$initial(),
    matrix(1, nrow = 2, ncol = 2, dimnames = list(c("p1", "p2"), c("b", "a")))
  )
})

test_that("pmcmc_nested_parameters names/populations", {
  parameters <- list(a = pmcmc_varied_parameter("a", c("p1", "p2"), 1:2),
                     b = pmcmc_varied_parameter("b", c("p1", "p2"), 3:4),
                     c = pmcmc_parameter("c", 5),
                     d = pmcmc_parameter("d", 6))
  proposal_fixed <- diag(2)
  proposal_varied <- diag(2) + 1
  p <- pmcmc_nested_parameters$new(parameters, proposal_varied, proposal_fixed)
  expect_equal(p$names(), c("c", "d", "a", "b"))
  expect_equal(p$populations(), c("p1", "p2"))

  # 1 varied 1 pop
  parameters <- list(a = pmcmc_varied_parameter("a", c("p1"), 1))
  proposal_varied <- diag(1)
  p <- pmcmc_nested_parameters$new(parameters, proposal_varied)
  expect_equal(p$names(), "a")
  expect_equal(p$populations(), "p1")
})

test_that("pmcmc_nested_parameters summary", {
  parameters <- list(a = pmcmc_varied_parameter("a", c("p1", "p2"), 1:2),
    b = pmcmc_varied_parameter("b", c("p1", "p2"), 3:4,
        discrete = c(TRUE, FALSE)),
    c = pmcmc_parameter("c", 5),
    d = pmcmc_parameter("d", 6))
  proposal_fixed <- diag(2)
  proposal_varied <- diag(2) + 1
  p <- pmcmc_nested_parameters$new(parameters, proposal_varied, proposal_fixed)
  expect_equal(p$summary("p1"),
    data.frame(name = c("c", "d", "a", "b"), min = -Inf, max = Inf,
        discrete = c(FALSE, FALSE, FALSE, TRUE),
        type = c("fixed", "fixed", "varied", "varied")))
  expect_equal(p$summary(), list(p1 = p$summary("p1"), p2 = p$summary("p2")))

  expect_error(p$summary("a"), "Expected 'population' in")

  # 1 varied 1 pop
  parameters <- list(a = pmcmc_varied_parameter("a", c("p1"), 1))
  proposal_varied <- diag(1)
  p <- pmcmc_nested_parameters$new(parameters, proposal_varied)
  expect_equal(p$summary(),
    list(p1 = data.frame(name = "a", min = -Inf, max = Inf, discrete = FALSE,
      type = "varied")))
})

test_that("pmcmc_nested_parameters prior", {
  parameters <- list(
    a = pmcmc_varied_parameter("a", c("p1", "p2"), 1:2,
        prior = list(function(x) 1, function(x) 2)),
    b = pmcmc_varied_parameter("b", c("p1", "p2"), 3:4,
        prior = list(function(x) 3, function(x) 4)),
    c = pmcmc_parameter("c", 5, prior = function(x) 5),
    d = pmcmc_parameter("d", 6, prior = function(x) 6))
  proposal_fixed <- diag(2)
  proposal_varied <- diag(2) + 1
  p <- pmcmc_nested_parameters$new(parameters, proposal_varied, proposal_fixed)

  init <- p$initial()
  expect_equal(p$prior(init), set_names(c(15, 17), c("p1", "p2")))

  parameters <- list(
    a = pmcmc_varied_parameter("a", c("p1", "p2"), 1:2,
        prior = list(dnorm, dexp)),
    b = pmcmc_varied_parameter("b", c("p1", "p2"), 3:4,
        prior = dlnorm),
    c = pmcmc_parameter("c", 5, prior = function(x) 5),
    d = pmcmc_parameter("d", 6, prior = function(x) 6))

  p <- pmcmc_nested_parameters$new(parameters, proposal_varied, proposal_fixed)

  init <- p$initial()
  expect_equal(p$prior(init),
    set_names(c(11 + dnorm(1) + dlnorm(3), 11 + dexp(2) + dlnorm(4)),
        c("p1", "p2")))

  # 1 varied 1 pop
  parameters <- list(a = pmcmc_varied_parameter("a", "p1", 1, prior = dnorm))
  proposal_varied <- diag(1)
  p <- pmcmc_nested_parameters$new(parameters, proposal_varied)
  init <- p$initial()
  expect_equal(p$prior(init), set_names(dnorm(1), "p1"))

  # 1 varied 1 fix
  parameters <- list(
    a = pmcmc_varied_parameter("a", c("p1", "p2"), 1, prior = dnorm),
    b = pmcmc_parameter("b", 1, prior = dexp))
  proposal_varied <- diag(1)
  proposal_fixed <- diag(1)
  p <- pmcmc_nested_parameters$new(parameters, proposal_varied, proposal_fixed)
  init <- p$initial()
  prior <- dnorm(1) + dexp(1)
  expect_equal(p$prior(init), set_names(rep(prior, 2), c("p1", "p2")))
})

test_that("pmcmc_nested_parameters model/transform", {
  parameters <- list(
    a = pmcmc_varied_parameter("a", c("p1", "p2"), 1:2,
        prior = list(function(x) 1, function(x) 2)),
    b = pmcmc_varied_parameter("b", c("p1", "p2"), 3:4,
        prior = list(function(x) 3, function(x) 4)),
    c = pmcmc_parameter("c", 5, prior = function(x) 5),
    d = pmcmc_parameter("d", 6, prior = function(x) 6))
  proposal_fixed <- diag(2)
  proposal_varied <- diag(2) + 1
  p <- pmcmc_nested_parameters$new(parameters, proposal_varied, proposal_fixed)

  expect_equal(p$model(p$initial()), apply(p$initial(), 1, as.list))

  p <- pmcmc_nested_parameters$new(parameters, proposal_varied, proposal_fixed,
    transform = function(x) as.list(log(x)))

  expect_equal(p$model(p$initial()),
                apply(p$initial(), 1, function(x) as.list(log(x))))
})

test_that("pmcmc_nested_parameters propose", {
  parameters <- list(
    a = pmcmc_varied_parameter("a", c("p1", "p2"), 1:2,
        prior = list(function(x) 1, function(x) 2)),
    b = pmcmc_varied_parameter("b", c("p1", "p2"), 3:4,
        prior = list(function(x) 3, function(x) 4)),
    c = pmcmc_parameter("c", 5, prior = function(x) 5),
    d = pmcmc_parameter("d", 6, prior = function(x) 6))
  proposal_fixed <- diag(2)
  proposal_varied <- diag(2) + 1
  p <- pmcmc_nested_parameters$new(parameters, proposal_varied, proposal_fixed)

  init <- p$initial()

  prop <- p$propose(init, type = "fixed")
  expect_true(identical(prop[, 3:4], init[, 3:4]))
  expect_false(identical(prop[, 1:2], init[, 1:2]))
  expect_equal(prop[1, 1:2], prop[2, 1:2])

  prop <- p$propose(init, type = "varied")
  expect_false(identical(prop[, 3:4], init[, 3:4]))
  expect_true(identical(prop[, 1:2], init[, 1:2]))
  expect_equal(prop[1, 1:2], prop[2, 1:2])

  prop <- p$propose(init, type = "both")
  expect_false(identical(prop[, 3:4], init[, 3:4]))
  expect_false(identical(prop[, 1:2], init[, 1:2]))
  expect_equal(prop[1, 1:2], prop[2, 1:2])

  # 1 varied 1 pop
  parameters <- list(a = pmcmc_varied_parameter("a", "p1", 1, prior = dnorm))
  proposal_varied <- diag(1)
  p <- pmcmc_nested_parameters$new(parameters, proposal_varied)
  init <- p$initial()

  prop <- p$propose(init, type = "varied")
  expect_true(inherits(prop, "matrix"))
  expect_identical(dimnames(prop), list("p1", "a"))

  expect_equal(p$propose(init, type = "fixed"), NULL)
  expect_equal(
    with(set.seed(1), p$propose(init, type = "both")),
    with(set.seed(1), p$propose(init, type = "varied"))
  )

  # 1 varied 1 fix
  parameters <- list(
    a = pmcmc_varied_parameter("a", "p1", 1, prior = dnorm),
    b = pmcmc_parameter("b", 1, prior = dexp))
  proposal_varied <- diag(1)
  proposal_fixed <- diag(1)
  p <- pmcmc_nested_parameters$new(parameters, proposal_varied, proposal_fixed)
  init <- p$initial()

  prop <- p$propose(init, type = "fixed")
  expect_true(identical(prop[, 2], init[, 2]))
  expect_false(identical(prop[, 1], init[, 1]))
})

# TODO
# test_that("pmcmc_nested_parameters fixed", {
#   p <- pmcmc_parameters$new(
#     parameters = list(
#       pmcmc_varied_parameter(name = "a", 2:4, min = 0:2, max = 4:6,
#         prior = list(dnorm, dexp, function(x) dnorm(x, 2))),
#       pmcmc_parameter(name = "b", 2, min = 1, max = 3, prior = dlnorm),
#       pmcmc_parameter(name = "c", 3),
#       pmcmc_varied_parameter(name = "d", 4)),
#     proposal = diag(4),
#     populations = c("Europe", "America", "Africa")
#   )
#   p_fix <- pmcmc_parameters$new(
#     parameters = list(
#       pmcmc_parameter(name = "b", 2, min = 1, max = 3, prior = dlnorm),
#       pmcmc_varied_parameter(name = "d", 4)),
#     proposal = diag(2),
#     populations = c("Europe", "America", "Africa"),
#     transform = function(p) {
#           base_transform(set_into(base, idx_vary, p))
#       }
#   )
#   fixedp <- p$fix(c(a = 1, c = 2))
#   expect_equal(fixedp, p_fix)

#   expect_equal(
#     fixedp$initial(),
#     matrix(c(2, 4), nrow = 3, ncol = 2, TRUE,
#     list(c("Europe", "America", "Africa"), c("b", "d")))
#   )

#   mod <- fixedp$model(fixedp$initial())
#   expect_true(is.list(mod) && length(mod) == 3)
#   expect_equal(names(mod), c("Europe", "America", "Africa"))
#   expect_equal(mod$Europe, list(a = 1, b = 2, c = 2, d = 4))
# })
