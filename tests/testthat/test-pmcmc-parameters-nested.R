
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

test_that("recycle", {
  expect_error(pmcmc_varied_parameter("p1", letters[1:3], 1:2),
               "Invalid length")
})

test_that("clean_parameters - error on empty and wrong type", {
  expect_error(clean_parameters(list()), "At least one")
  expect_error(clean_parameters(list("a")), "Expected all")
})

test_that("clean_parameters - error on missing proposal", {
  expect_error(clean_parameters(
    list(pmcmc_varied_parameter("a", "a", 2),
         pmcmc_varied_parameter("a", "a", 3)), NULL, NULL),
    "'proposal_varied'")
  expect_error(clean_parameters(
    list(pmcmc_parameter("a", 2),
         pmcmc_parameter("a", 3)), NULL, NULL),
    "'proposal_fixed'")
})

test_that("clean_parameters - error on duplicates", {
  expect_error(clean_parameters(
    list(pmcmc_varied_parameter("a", "a", 2),
         pmcmc_varied_parameter("a", "a", 3)), diag(2), diag(2)),
    "Duplicate")
  expect_error(clean_parameters(
    list(pmcmc_varied_parameter("a", "a", 2),
         pmcmc_parameter("a", 3)), diag(2), diag(2)),
    "Duplicate")
})

test_that("clean_parameters - error on wrong names", {
  expect_error(clean_parameters(
    list(a = pmcmc_varied_parameter("a", "a", 2),
         c = pmcmc_parameter("b", 3)), diag(2), diag(2)),
    "Fixed parameters are named")
  expect_error(clean_parameters(
    list(c = pmcmc_varied_parameter("a", "a", 2),
         b = pmcmc_parameter("b", 3)), diag(2), diag(2)),
    "Varied parameters are named")
  expect_error(clean_parameters(
    list(pmcmc_varied_parameter("a", c("p1", "p2"), 2),
         pmcmc_varied_parameter("b", c("p2", "p1"), 2),
         pmcmc_parameter("c", 3)), diag(2), diag(2)),
    "Populations and ordering")
  expect_error(clean_parameters(
    list(pmcmc_varied_parameter("a", c("p1", "p2"), 2),
         pmcmc_varied_parameter("b", c("p2"), 2),
         pmcmc_parameter("c", 3)), diag(2), diag(2)),
    "Populations and ordering")
})

test_that("clean_parameters - silent", {
  p <- clean_parameters(
    list(a = pmcmc_varied_parameter("a", c("p1", "p2"), 2),
         b = pmcmc_parameter("b", 3),
         c = pmcmc_parameter("c", 3),
         d = pmcmc_varied_parameter("d", c("p1", "p2"), 2)),
    diag(2), diag(2))
  expect_equal(p$names,
               list(fixed = c("b", "c"), varied = c("a", "d"),
                    original = letters[1:4]))
})

test_that("clean_parameters - silent without fixed/varied", {
  expect_equal(clean_parameters(
    list(a = pmcmc_varied_parameter("a", c("p1", "p2"), 2),
         b = pmcmc_varied_parameter("b", c("p1", "p2"), 2)),
    diag(2), NULL)$fixed_parameters, NULL)
  expect_equal(clean_parameters(
    list(a = pmcmc_parameter("a", 2),
         b = pmcmc_parameter("b", 2)), NULL, diag(1))$varied_parameters, NULL)
})

test_that("clean proposals - error on wrong names", {
  params <- list(fixed = letters[1:2], varied = letters[3:4])
  pops <- letters[5:6]

  expect_error(
    clean_proposals(array(diag(2), c(2, 2, 2),
                          dimnames = list(letters[1:2], letters[3:4],
                                          letters[5:6])),
                    diag(2), pops, params),
    "proposal_varied")

  expect_error(
    clean_proposals(
      array(diag(2), c(2, 2, 2)),
      matrix(diag(2), 2, 2, dimnames = list(letters[5:6], letters[5:6])),
      pops, params),
    "proposal_fixed")

  expect_error(
    clean_proposals(array(diag(2), c(2, 2, 2),
                          dimnames = list(NULL, NULL, letters[3:4])),
                    diag(2), pops, params),
    "Expected 3rd")
})

test_that("clean proposals - error on misspecified array", {
  params <- list(fixed = letters[1:2], varied = letters[3:4])
  pops <- letters[5:6]

  expect_error(
    clean_proposals(array(diag(2), c(1, 2, 2),
                          dimnames = list(NULL, NULL, letters[3:4])),
                    diag(2), pops, params),
    "Expected proposal")
})

test_that("clean proposals - silent", {
  params <- list(fixed = letters[1:2], varied = letters[3:4])
  pops <- letters[5:6]

  expect_equal(
    clean_proposals(diag(2), diag(2), pops, params)$
      varied,
    array(diag(2), dim = c(2, 2, 2),
          dimnames = list(letters[3:4], letters[3:4], letters[5:6])))

  expect_equal(
    clean_proposals(diag(2), diag(2), pops, params)$
      fixed,
    matrix(diag(2), 2, 2, dimnames = list(letters[1:2], letters[1:2])))
})

test_that("clean proposals - silent edge cases", {
  params <- list(fixed = letters[1:2], varied = letters[3:4])
  pops <- letters[5:6]

  # 1 varied 1 fixed
  expect_equal(
    clean_proposals(diag(1), diag(1), pops, list(fixed = "a", varied = "b"))$
      varied,
    array(diag(1), dim = c(1, 1, 2),
          dimnames = list("b", "b", pops)))

  # 1 varied 0 fixed 1 pop
  expect_equal(
    clean_proposals(diag(1), NULL, "c", list(varied = "b"))$
      varied,
    array(diag(1), dim = c(1, 1, 1),
          dimnames = list("b", "b", "c")))

  # 0 varied 1 fixed
  expect_equal(
    clean_proposals(NULL, diag(1), NULL, list(fixed = "b"))$
      fixed,
    matrix(diag(1), dimnames = list("b", "b")))
})

test_that("make proposals - varied and fixed", {
  params <- list(fixed = letters[1:2], varied = letters[3:4])
  pops <- letters[5:6]
  proposal_varied <- diag(2)
  proposal_fixed <- diag(2)
  kernels <- clean_proposals(proposal_varied, proposal_fixed, pops, params)

  props <- make_proposals(kernels$varied, kernels$fixed,
                          list(original = letters[1:4]))
  expect_equal(props$kernel_varied$e,
               matrix(c(rep(0, 10), 1, rep(0, 4), 1), 4, 4,
                      dimnames = rep(list(letters[1:4]), 2)))
  expect_equal(props$kernel_fixed, kernels$fixed)
  expect_true(inherits(props$proposal_varied$e, "function"))
  expect_true(inherits(props$proposal_varied$f, "function"))
})

test_that("make proposals - 1 varied, 1 fixed", {
  params <- list(fixed = "a", varied = "b")
  pops <- "p"
  proposal_varied <- diag(1)
  proposal_fixed <- diag(1)
  kernels <- clean_proposals(proposal_varied, proposal_fixed, pops, params)

  props <- make_proposals(kernels$varied, kernels$fixed,
                          list(original = letters[1:2]))
  expect_equal(props$kernel_varied$p,
               matrix(c(rep(0, 3), 1), 2, 2,
                      dimnames = rep(list(letters[1:2]), 2)))
  expect_equal(props$kernel_fixed, kernels$fixed)
  expect_true(inherits(props$proposal_varied$p, "function"))
  expect_true(inherits(props$proposal_varied$p, "function"))
})

test_that("make proposals - varied only", {
  params <- list(varied = letters[3:4])
  pops <- letters[5:6]
  proposal_varied <- diag(2)
  kernels <- clean_proposals(proposal_varied, NULL, pops, params)

  props <- make_proposals(kernels$varied, kernels$fixed,
                          list(original = letters[3:4]))
  expect_equal(props$kernel_varied$e, kernels$varied[, , 1])
  expect_equal(props$kernel_fixed, kernels$fixed)
  expect_true(inherits(props$proposal_varied$e, "function"))
  expect_true(inherits(props$proposal_varied$f, "function"))

  props <- make_proposals(kernels$varied, NULL,
                          list(original = letters[1:4]))
  expect_equal(props$proposal_fixed, NULL)
})

test_that("make proposals - fixed only", {
  params <- list(fixed = letters[3:4])
  pops <- letters[5:6]
  proposal_fixed <- diag(2)
  kernels <- clean_proposals(NULL, proposal_fixed, NULL, params)

  props <- make_proposals(kernels$varied, kernels$fixed, letters[3:4])
  expect_equal(props$kernel_varied$e, kernels$varied[, , 1])
  expect_equal(props$kernel_fixed, kernels$fixed)
  expect_true(inherits(props$proposal_fixed, "function"))
  expect_equal(props$proposal_varied, NULL)
})

test_that("construct pmcmc_parameters_nested - silent", {
  parameters <- list(a = pmcmc_varied_parameter("a", c("p1", "p2"), 2),
                     b = pmcmc_varied_parameter("b", c("p1", "p2"), 2),
                     c = pmcmc_parameter("c", 3),
                     d = pmcmc_parameter("d", 4))
  proposal_fixed <- diag(2)
  proposal_varied <- diag(2) + 1

  expect_silent(pmcmc_parameters_nested$new(parameters, proposal_varied,
                                            proposal_fixed))
  expect_silent(pmcmc_parameters_nested$new(parameters, proposal_varied,
                                            proposal_fixed,
                                            transform = as.list))
})

test_that("construct pmcmc_parameters_nested - 1 varied 1 fixed 1 pop", {
  parameters <- list(a = pmcmc_varied_parameter("a", c("p2"), 2),
                     d = pmcmc_parameter("d", 4))
  proposal_fixed <- diag(1)
  proposal_varied <- diag(1) + 1
  expect_silent(pmcmc_parameters_nested$new(parameters, proposal_varied,
                                            proposal_fixed))
})

test_that("construct pmcmc_parameters_nested - 1 varied 0 fixed", {
  parameters <- list(a = pmcmc_varied_parameter("a", c("p1", "p2"), 2))
  proposal_varied <- diag(1) + 1
  expect_silent(pmcmc_parameters_nested$new(parameters, proposal_varied))
})

test_that("construct pmcmc_parameters_nested - 0 varied 1 fixed", {
  parameters <- list(a = pmcmc_parameter("a", 2))
  proposal_fixed <- diag(1)
  expect_error(
    pmcmc_parameters_nested$new(parameters, NULL, proposal_fixed), "Either")
  expect_silent(
    pmcmc_parameters_nested$new(parameters, NULL, proposal_fixed, "a"))
})

test_that("pmcmc_parameters_nested initial", {
  parameters <- list(a = pmcmc_varied_parameter("a", c("p1", "p2"), 3:4),
                     b = pmcmc_varied_parameter("b", c("p1", "p2"), 1:2),
                     c = pmcmc_parameter("c", 5),
                     d = pmcmc_parameter("d", 6))
  proposal_fixed <- diag(2)
  proposal_varied <- diag(2) + 1
  p <- pmcmc_parameters_nested$new(parameters, proposal_varied, proposal_fixed)
  expect_equal(
    p$initial(),
    matrix(c(3:4, 1:2, 5, 5, 6, 6), nrow = 2, ncol = 4,
           dimnames = list(c("p1", "p2"), letters[1:4]))
  )
})

test_that("pmcmc_parameters_nested initial - varied only", {
  parameters <- list(
    a = pmcmc_varied_parameter("a", c("p1", "p2"), 1:2),
    b = pmcmc_varied_parameter("b", c("p1", "p2"), 3:4)
  )
  proposal_varied <- diag(2)
  p <- pmcmc_parameters_nested$new(parameters, proposal_varied)
  expect_equal(
    p$initial(),
    matrix(1:4, 2, 2, dimnames = list(c("p1", "p2"), letters[1:2]))
  )

  parameters <- list(
    a = pmcmc_varied_parameter("a", c("p1"), 1)
  )
  proposal_varied <- diag(1)
  p <- pmcmc_parameters_nested$new(parameters, proposal_varied)
  expect_equal(
    p$initial(),
    matrix(1, nrow = 1, ncol = 1, dimnames = list(c("p1"), "a"))
  )
})

test_that("pmcmc_parameters_nested initial - fixed only", {
  parameters <- list(
    a = pmcmc_parameter("a", 1),
    b = pmcmc_parameter("b", 3)
  )
  proposal_fixed <- diag(2)
  p <- pmcmc_parameters_nested$new(parameters, NULL, proposal_fixed,
                                   c("p1", "p2"))
  expect_equal(
    p$initial(),
    matrix(c(1, 1, 3, 3), 2, 2, dimnames = list(c("p1", "p2"), letters[1:2]))
  )
})

test_that("pmcmc_parameters_nested initial - 1 varied 1  fix", {
  parameters <- list(
    a = pmcmc_varied_parameter("a", c("p1", "p2"), 1, prior = dnorm),
    b = pmcmc_parameter("b", 2, prior = dexp))
  proposal_varied <- diag(1)
  proposal_fixed <- diag(1)
  p <- pmcmc_parameters_nested$new(parameters, proposal_varied, proposal_fixed)
  expect_equal(
    p$initial(),
    matrix(c(1, 1, 2, 2), nrow = 2, ncol = 2,
           dimnames = list(c("p1", "p2"), c("a", "b")))
  )
})

test_that("pmcmc_parameters_nested names/populations", {
  parameters <- list(a = pmcmc_varied_parameter("a", c("p1", "p2"), 1:2),
                     b = pmcmc_varied_parameter("b", c("p1", "p2"), 3:4),
                     c = pmcmc_parameter("c", 5),
                     d = pmcmc_parameter("d", 6))
  proposal_fixed <- diag(2)
  proposal_varied <- diag(2) + 1
  p <- pmcmc_parameters_nested$new(parameters, proposal_varied, proposal_fixed)
  expect_equal(p$names(), letters[1:4])
  expect_equal(p$populations(), c("p1", "p2"))
})

test_that("pmcmc_parameters_nested summary", {
  parameters <- list(
    a = pmcmc_varied_parameter("a", c("p1", "p2"), 1:2),
    b = pmcmc_parameter("b", 5),
    c = pmcmc_varied_parameter("c", c("p1", "p2"), 3:4,
                               discrete = c(TRUE, FALSE)),
    d = pmcmc_parameter("d", 6))
  proposal_fixed <- diag(2)
  proposal_varied <- diag(2) + 1
  p <- pmcmc_parameters_nested$new(parameters, proposal_varied, proposal_fixed)

  expect_equal(p$summary("p1"),
               data.frame(name = letters[1:4], min = -Inf, max = Inf,
                          discrete = c(FALSE, FALSE, TRUE, FALSE),
                          type = c("varied", "fixed", "varied", "fixed")))
  expect_equal(p$summary(), list(p1 = p$summary("p1"), p2 = p$summary("p2")))

  expect_error(p$summary("a"), "Expected 'population' in")
})

test_that("pmcmc_parameters_nested summary - 1 varied 1 pop 0 fix", {
  parameters <- list(a = pmcmc_varied_parameter("a", c("p1"), 1))
  proposal_varied <- diag(1)
  p <- pmcmc_parameters_nested$new(parameters, proposal_varied)
  expect_equal(p$summary(),
               list(p1 = data.frame(name = "a", min = -Inf, max = Inf,
                                    discrete = FALSE, type = "varied")))
})

test_that("pmcmc_parameters_nested summary - 0 varied 1 pop 1 fix", {
  parameters <- list(a = pmcmc_parameter("a", 1))
  proposal_fixed <- diag(1)
  p <- pmcmc_parameters_nested$new(parameters, NULL, proposal_fixed, "p1")
  expect_equal(p$summary(),
               list(p1 = data.frame(name = "a", min = -Inf, max = Inf,
                                    discrete = FALSE, type = "fixed")))
})

test_that("pmcmc_parameters_nested prior", {
  parameters <- list(
    a = pmcmc_varied_parameter("a", c("p1", "p2"), 1:2,
                               prior = list(function(x) 1, function(x) 2)),
    b = pmcmc_varied_parameter("b", c("p1", "p2"), 3:4,
                               prior = list(function(x) 3, function(x) 4)),
    c = pmcmc_parameter("c", 5, prior = function(x) 5),
    d = pmcmc_parameter("d", 6, prior = function(x) 6))
  proposal_fixed <- diag(2)
  proposal_varied <- diag(2) + 1
  p <- pmcmc_parameters_nested$new(parameters, proposal_varied, proposal_fixed)

  init <- p$initial()
  expect_equal(p$prior(init), set_names(c(15, 17), c("p1", "p2")))

  parameters <- list(
    a = pmcmc_varied_parameter("a", c("p1", "p2"), 1:2,
                               prior = list(dnorm, dexp)),
    b = pmcmc_varied_parameter("b", c("p1", "p2"), 3:4,
                               prior = dlnorm),
    c = pmcmc_parameter("c", 5, prior = function(x) 5),
    d = pmcmc_parameter("d", 6, prior = function(x) 6))

  p <- pmcmc_parameters_nested$new(parameters, proposal_varied, proposal_fixed)

  init <- p$initial()
  expect_equal(p$prior(init),
               set_names(c(11 + dnorm(1) + dlnorm(3), 11 + dexp(2) + dlnorm(4)),
                         c("p1", "p2")))
})

test_that("pmcmc_parameters_nested prior - 1 varied 1 pop 0 fix", {
  parameters <- list(a = pmcmc_varied_parameter("a", "p1", 1, prior = dnorm))
  proposal_varied <- diag(1)
  p <- pmcmc_parameters_nested$new(parameters, proposal_varied)
  init <- p$initial()
  expect_equal(p$prior(init), set_names(dnorm(1), "p1"))
})

test_that("pmcmc_parameters_nested prior - 1 varied 1 fix", {
  parameters <- list(
    a = pmcmc_varied_parameter("a", c("p1", "p2"), 1, prior = dnorm),
    b = pmcmc_parameter("b", 1, prior = dexp))
  proposal_varied <- diag(1)
  proposal_fixed <- diag(1)
  p <- pmcmc_parameters_nested$new(parameters, proposal_varied, proposal_fixed)
  init <- p$initial()
  prior <- dnorm(1) + dexp(1)
  expect_equal(p$prior(init), set_names(rep(prior, 2), c("p1", "p2")))
})

test_that("pmcmc_parameters_nested prior - 0 varied 1 fix", {
  parameters <- list(
    b = pmcmc_parameter("b", 1, prior = dexp))
  proposal_fixed <- diag(1)
  p <- pmcmc_parameters_nested$new(parameters, NULL, proposal_fixed,
                                   c("p1", "p2"))
  init <- p$initial()
  expect_equal(p$prior(init), set_names(rep(dexp(1), 2), c("p1", "p2")))
})

test_that("pmcmc_parameters_nested model/transform", {
  parameters <- list(
    a = pmcmc_varied_parameter("a", c("p1", "p2"), 1:2,
                               prior = list(function(x) 1, function(x) 2)),
    b = pmcmc_varied_parameter("b", c("p1", "p2"), 3:4,
                               prior = list(function(x) 3, function(x) 4)),
    c = pmcmc_parameter("c", 5, prior = function(x) 5),
    d = pmcmc_parameter("d", 6, prior = function(x) 6))
  proposal_fixed <- diag(2)
  proposal_varied <- diag(2) + 1
  p <- pmcmc_parameters_nested$new(parameters, proposal_varied, proposal_fixed)

  expect_equal(p$model(p$initial()), apply(p$initial(), 1, as.list))

  p <- pmcmc_parameters_nested$new(parameters, proposal_varied, proposal_fixed,
                                   transform = function(x)
                                     apply(x, 1, function(y) as.list(log(y))))

  expect_equal(p$model(p$initial()),
               apply(p$initial(), 1, function(x) as.list(log(x))))
})

test_that("pmcmc_parameters_nested propose", {
  parameters <- list(
    a = pmcmc_varied_parameter("a", c("p1", "p2"), 1:2,
                               prior = list(function(x) 1, function(x) 2)),
    b = pmcmc_parameter("b", 5, prior = function(x) 5),
    c = pmcmc_varied_parameter("c", c("p1", "p2"), 3:4,
                               prior = list(function(x) 3, function(x) 4)),
    d = pmcmc_parameter("d", 6, prior = function(x) 6))
  proposal_fixed <- diag(2)
  proposal_varied <- diag(2) + 1
  p <- pmcmc_parameters_nested$new(parameters, proposal_varied, proposal_fixed)

  init <- p$initial()

  prop <- p$propose(init, type = "fixed")
  expect_identical(colnames(prop), letters[1:4])
  expect_false(identical(prop[, c(2, 4)], init[, c(2, 4)]))
  expect_true(identical(prop[, c(1, 3)], init[, c(1, 3)]))
  expect_equal(prop[1, c(2, 4)], prop[2, c(2, 4)])

  prop <- p$propose(init, type = "varied")
  expect_identical(colnames(prop), letters[1:4])
  expect_true(identical(prop[, c(2, 4)], init[, c(2, 4)]))
  expect_false(identical(prop[, c(1, 3)], init[, c(1, 3)]))
  expect_equal(prop[1, c(2, 4)], prop[2, c(2, 4)])

  prop <- p$propose(init, type = "both")
  expect_identical(colnames(prop), letters[1:4])
  expect_false(identical(prop[, c(2, 4)], init[, c(2, 4)]))
  expect_false(identical(prop[, c(1, 3)], init[, c(1, 3)]))
  expect_equal(prop[1, c(2, 4)], prop[2, c(2, 4)])
})

test_that("pmcmc_parameters_nested propose - 1 varied 1 pop", {
  parameters <- list(a = pmcmc_varied_parameter("a", "p1", 1, prior = dnorm))
  proposal_varied <- diag(1)
  p <- pmcmc_parameters_nested$new(parameters, proposal_varied)
  init <- p$initial()

  prop <- p$propose(init, type = "varied")
  expect_true(inherits(prop, "matrix"))
  expect_identical(dimnames(prop), list("p1", "a"))

  expect_equal(p$propose(init, type = "fixed"), NULL)
  expect_equal(
    with(set.seed(1), p$propose(init, type = "both")),
    with(set.seed(1), p$propose(init, type = "varied"))
  )
})

test_that("pmcmc_parameters_nested propose - 1 varied 1 fix", {
  parameters <- list(
    a = pmcmc_varied_parameter("a", "p1", 1, prior = dnorm),
    b = pmcmc_parameter("b", 1, prior = dexp))
  proposal_varied <- diag(1)
  proposal_fixed <- diag(1)
  p <- pmcmc_parameters_nested$new(parameters, proposal_varied, proposal_fixed)
  init <- p$initial()

  prop <- p$propose(init, type = "fixed")
  expect_true(identical(prop[, 1], init[, 1]))
  expect_false(identical(prop[, 2], init[, 2]))
})


test_that("pmcmc_parameters_nested fixed", {
  parameters <- list(
    a = pmcmc_varied_parameter("a", c("p1", "p2"), 1:2,
                               prior = list(function(x) 1, function(x) 2)),
    b = pmcmc_parameter("b", 1,
                        prior = function(x) 1),
    c = pmcmc_varied_parameter("c", c("p1", "p2"), 1:2,
                               prior = list(function(x) 1, function(x) 2)),
    d = pmcmc_parameter("d", 5, prior = function(x) 5),
    e = pmcmc_varied_parameter("e", c("p1", "p2"), 3:4,
                               prior = list(function(x) 3, function(x) 4)),
    f = pmcmc_parameter("f", 6, prior = function(x) 6))
  proposal_fixed <- diag(3)
  proposal_varied <- diag(3) + 1
  p <- pmcmc_parameters_nested$new(parameters, proposal_varied, proposal_fixed)

  fixed <- matrix(c(1, 1, 2, 2, 3, 4), 2, 3,
                  dimnames = list(c("p1", "p2"), letters[1:3]))
  fix_p <- expect_silent(p$fix(fixed))

  fixed_p <- pmcmc_parameters_nested$new(
    list(
      d = pmcmc_parameter("d", 5, prior = function(x) 5),
      e = pmcmc_varied_parameter("e", c("p1", "p2"), 3:4,
                                 prior = list(function(x) 3, function(x) 4)),
      f = pmcmc_parameter("f", 6, prior = function(x) 6)),
    diag(1) + 1, diag(2),
    transform = function(p) {
      base[, idx_vary] <- p
      base_transform(base)
    })

  expect_equal(fix_p, fixed_p)

  init <- fix_p$initial()
  expect_identical(fix_p$model(init),
                   list(p1 = list(a = 1, b = 2, c = 3, d = 5, e = 3, f = 6),
                        p2 = list(a = 1, b = 2, c = 4, d = 5, e = 4, f = 6))
  )
})

test_that("pmcmc_parameters_nested fixed - 1 fix 1 var 1 pop - fix vary", {
  parameters <- list(
    a = pmcmc_varied_parameter("a", c("p1"), 1,
                               prior = function(x) 2),
    b = pmcmc_parameter("b", 1,
                        prior = function(x) 1))
  proposal_fixed <- diag(1)
  proposal_varied <- diag(1) + 1
  p <- pmcmc_parameters_nested$new(parameters, proposal_varied, proposal_fixed)

  fixed <- matrix(1, dimnames = list("p1", "a"))
  fix_p <- expect_silent(p$fix(fixed))

  fixed_p <- pmcmc_parameters_nested$new(
    list(b = pmcmc_parameter("b", 1, prior = function(x) 1)),
    NULL, diag(1), "p1",
    transform = function(p) {
      base[, idx_vary] <- p
      base_transform(base)
    })

  expect_equal(fix_p, fixed_p)

  init <- fix_p$initial()
  expect_identical(fix_p$model(init),
                   list(p1 = list(a = 1, b = 1))
  )
})

test_that("pmcmc_parameters_nested fixed - 1 fix 1 var 1 pop - fix fixed", {
  parameters <- list(
    a = pmcmc_varied_parameter("a", c("p1"), 1,
                               prior = function(x) 2),
    b = pmcmc_parameter("b", 1,
                        prior = function(x) 1))
  proposal_fixed <- diag(1)
  proposal_varied <- diag(1) + 1
  p <- pmcmc_parameters_nested$new(parameters, proposal_varied, proposal_fixed)

  fixed <- matrix(3, dimnames = list("p1", "b"))
  fix_p <- expect_silent(p$fix(fixed))

  fixed_p <- pmcmc_parameters_nested$new(
    list(a = pmcmc_varied_parameter("a", c("p1"), 1, prior = function(x) 2)),
    diag(1) + 1, NULL,
    transform = function(p) {
      base[, idx_vary] <- p
      base_transform(base)
    })
  r6_private(r6_private(fix_p)$varied_parameters$p1)$proposal_kernel
  r6_private(r6_private(fixed_p)$varied_parameters$p1)$proposal_kernel
  expect_equal(fix_p, fixed_p)

  init <- fix_p$initial()
  expect_identical(fix_p$model(init),
                   list(p1 = list(a = 1, b = 3))
  )
})
