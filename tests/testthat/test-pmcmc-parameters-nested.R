test_that("Can construct a varied parameter", {
  p <- pmcmc_varied_parameter("p1", letters[1:2], 1:2)
  expect_s3_class(p, "pmcmc_varied_parameter")
  expect_equal(p$a$name, "p1")
  expect_equal(p$b$name, "p1")
  expect_equal(p$a$initial, 1)
  expect_equal(p$b$initial, 2)
  expect_equal(p$a$min, -Inf)
  expect_equal(p$b$min, -Inf)
  expect_equal(p$a$max, Inf)
  expect_equal(p$b$max, Inf)
  expect_equal(p$a$discrete, FALSE)
  expect_equal(p$b$discrete, FALSE)
  expect_equal(p$a$prior(1), 0)
  expect_equal(p$b$prior(1), 0)
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


test_that("construct parameters - error on empty and wrong type", {
  expect_error(
    pmcmc_parameters_nested$new(list()),
    "At least one parameter is required")
  expect_error(
    pmcmc_parameters_nested$new(list(TRUE)),
    paste("Elements of 'parameters' must be in",
          "'{pmcmc_parameter, pmcmc_varied_parameter}'"),
    fixed = TRUE)
})


test_that("construct parameters - error on missing proposal", {
  pars_varied <- list(pmcmc_varied_parameter("a", "x", 2),
                      pmcmc_varied_parameter("b", "x", 3))
  pars_fixed <- list(pmcmc_parameter("c", 2),
                     pmcmc_parameter("d", 3))
  expect_error(
    pmcmc_parameters_nested$new(pars_varied),
    "'proposal_varied' not supplied for varied parameters")
  expect_error(
    pmcmc_parameters_nested$new(pars_fixed, populations = "x"),
    "'proposal_fixed' not supplied for fixed parameters")
})


test_that("construct parameters - error on unwanted proposal", {
  pars_varied <- list(pmcmc_varied_parameter("a", "x", 2),
                      pmcmc_varied_parameter("b", "x", 3))
  pars_fixed <- list(pmcmc_parameter("c", 2),
                     pmcmc_parameter("d", 3))
  expect_error(
    pmcmc_parameters_nested$new(pars_varied, diag(2), diag(2)),
    "'proposal_fixed' supplied, but no fixed parameters")
  expect_error(
    pmcmc_parameters_nested$new(pars_fixed, diag(2), diag(2), "x"),
    "'proposal_varied' supplied, but no varied parameters")
})


test_that("require explicit populations if no varied parameters", {
  pars_fixed <- list(pmcmc_parameter("c", 2),
                     pmcmc_parameter("d", 3))
  expect_error(
    pmcmc_parameters_nested$new(pars_fixed, NULL, diag(2)),
    paste("Either varied parameters must be included in 'parameters' or",
          "'populations' must be non-NULL"))
  p <- pmcmc_parameters_nested$new(pars_fixed, NULL, diag(2), c("x", "y"))
  expect_equal(p$populations(), c("x", "y"))
})


test_that("if explicit population provided, must match varied", {
  pars <- list(
    pmcmc_varied_parameter("a", "x", 2),
    pmcmc_varied_parameter("b", "x", 3),
    pmcmc_parameter("c", 2),
    pmcmc_parameter("d", 3))
  expect_error(
    pmcmc_parameters_nested$new(pars, diag(2), diag(2), "y"),
    "'population' does not match varied parameters")
})


test_that("construct parameters - error on duplicates", {
  pars_varied <- list(pmcmc_varied_parameter("a", "a", 2),
                      pmcmc_varied_parameter("a", "a", 3))
  pars_fixed <- list(pmcmc_parameter("a", 2),
                     pmcmc_parameter("a", 3))

  expect_error(
    pmcmc_parameters_nested$new(pars_varied, diag(2)),
    "Duplicate parameter names: 'a'")
  expect_error(
    pmcmc_parameters_nested$new(pars_fixed, NULL, diag(2)),
    "Duplicate parameter names: 'a'")
})


test_that("construct parameters - error on wrong names", {
  expect_error(
    pmcmc_parameters_nested$new(
      list(a = pmcmc_varied_parameter("a", "a", 2),
           c = pmcmc_parameter("b", 3)), diag(1), diag(1)),
    "Fixed parameters are named, but the names do not match parameters")
  expect_error(pmcmc_parameters_nested$new(
    list(c = pmcmc_varied_parameter("a", "a", 2),
         b = pmcmc_parameter("b", 3)), diag(1), diag(1)),
    "Varied parameters are named, but the names do not match parameters")
  expect_error(pmcmc_parameters_nested$new(
    list(pmcmc_varied_parameter("a", c("p1", "p2"), 2),
         pmcmc_varied_parameter("b", c("p2", "p1"), 2),
         pmcmc_parameter("c", 3)), diag(2), diag(1)),
    "Populations and ordering of varied parameters must be identical")
  expect_error(pmcmc_parameters_nested$new(
    list(pmcmc_varied_parameter("a", c("p1", "p2"), 2),
         pmcmc_varied_parameter("b", c("p2"), 2),
         pmcmc_parameter("c", 3)), diag(2), diag(1)),
    "Populations and ordering of varied parameters must be identical")
})


test_that("clean proposals - error on wrong names", {
  pars <- list(
    pmcmc_varied_parameter("a", c("x", "y"), 1:2),
    pmcmc_varied_parameter("b", c("x", "y"), 3:4),
    pmcmc_parameter("c", 5),
    pmcmc_parameter("d", 6))

  proposal_varied <- array(
    diag(2), c(2, 2, 2),
    dimnames = list(c("a", "b"), c("c", "d"), c("x", "y")))
  expect_error(
    pmcmc_parameters_nested$new(pars, proposal_varied, diag(2)),
    "Expected names of dimension 2 of 'proposal_varied' to match parameters")

  proposal_fixed <- matrix(diag(2), 2,
                           dimnames = list(c("a", "b"), c("c", "d")))
  expect_error(
    pmcmc_parameters_nested$new(pars, diag(2), proposal_fixed),
    "Expected names of dimension 1 of 'proposal_fixed' to match parameters")
})


test_that("clean proposals - error on misspecified array", {
  pars <- list(
    pmcmc_varied_parameter("a", c("x", "y"), 1:2),
    pmcmc_varied_parameter("b", c("x", "y"), 3:4),
    pmcmc_parameter("c", 5),
    pmcmc_parameter("d", 6))

  expect_error(
    pmcmc_parameters_nested$new(pars, diag(3), diag(2)),
    "Expected 'proposal_varied' to be array with dimensions 2 x 2")
  expect_error(
    pmcmc_parameters_nested$new(pars, array(1, c(2, 2, 3)), diag(2)),
    "Expected 'proposal_varied' to be array with dimensions 2 x 2 x 2")
  expect_error(
    pmcmc_parameters_nested$new(pars, diag(2), diag(3)),
    "Expected 'proposal_fixed' to be array with dimensions 2 x 2")
})


test_that("construct pmcmc_parameters_nested; contruction and basic use", {
  parameters <- list(a = pmcmc_varied_parameter("a", c("p1", "p2"), 1:2),
                     b = pmcmc_varied_parameter("b", c("p1", "p2"), 3:4),
                     c = pmcmc_parameter("c", 5),
                     d = pmcmc_parameter("d", 6))
  proposal_fixed <- diag(2)
  proposal_varied <- diag(2) + 1

  res <- pmcmc_parameters_nested$new(parameters, proposal_varied,
                                     proposal_fixed)
  expect_s3_class(res, "pmcmc_parameters_nested")
  expect_equal(
    res$initial(),
    cbind(p1 = c(a = 1, b = 3, c = 5, d = 6),
          p2 = c(a = 2, b = 4, c = 5, d = 6)))
  expect_equal(
    res$model(res$initial()),
    unname(apply(res$initial(), 2, as.list)))

  expect_equal(res$names(), c("a", "b", "c", "d"))
  expect_equal(res$names("fixed"), c("c", "d"))
  expect_equal(res$names("varied"), c("a", "b"))
  expect_equal(res$populations(), c("p1", "p2"))

  expect_equal(
    res$summary(),
    data_frame(name = rep(letters[1:4], 2),
               min = -Inf, max = Inf, discrete = FALSE,
               type = rep(c("varied", "fixed"), each = 2),
               population = rep(c("p1", "p2"), each = 4)))
})


test_that("pmcmc_parameters_nested initial - varied only", {
  parameters <- list(
    a = pmcmc_varied_parameter("a", c("p1", "p2"), 1:2),
    b = pmcmc_varied_parameter("b", c("p1", "p2"), 3:4))
  proposal_varied <- diag(2)
  p <- pmcmc_parameters_nested$new(parameters, proposal_varied)
  expect_equal(
    p$initial(),
    cbind(p1 = c(a = 1, b = 3), p2 = c(a = 2, b = 4)))

  expect_equal(p$names("fixed"), NULL)
  expect_equal(p$names("varied"), c("a", "b"))

  parameters <- list(
    a = pmcmc_varied_parameter("a", c("p1", "p2"), 1:2))
  proposal_varied <- diag(1)
  p <- pmcmc_parameters_nested$new(parameters, proposal_varied)
  expect_equal(
    p$initial(),
    cbind(p1 = c(a = 1), p2 = c(a = 2)))

  expect_equal(p$names("fixed"), NULL)
  expect_equal(p$names("varied"), "a")
})


test_that("pmcmc_parameters_nested initial - fixed only", {
  parameters <- list(
    a = pmcmc_parameter("a", 1),
    b = pmcmc_parameter("b", 3))
  proposal_fixed <- diag(2)
  p <- pmcmc_parameters_nested$new(parameters, NULL, proposal_fixed,
                                   c("p1", "p2"))
  expect_equal(
    p$initial(),
    cbind(p1 = c(a = 1, b = 3), p2 = c(a = 1, b = 3)))

  expect_equal(p$names("varied"), NULL)
  expect_equal(p$names("fixed"), c("a", "b"))
})


test_that("pmcmc_parameters_nested initial - 1 varied 1  fix", {
  parameters <- list(
    a = pmcmc_varied_parameter("a", c("p1", "p2"), 1:2, prior = dnorm),
    b = pmcmc_parameter("b", 3, prior = dexp))
  proposal_varied <- diag(1)
  proposal_fixed <- diag(1)
  p <- pmcmc_parameters_nested$new(parameters, proposal_varied, proposal_fixed)
  expect_equal(
    p$initial(),
    cbind(p1 = c(a = 1, b = 3), p2 = c(a = 2, b = 3)))

  expect_equal(p$names("varied"), "a")
  expect_equal(p$names("fixed"), "b")
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


test_that("pmcmc_parameters_nested model/transform", {
  parameters <- list(
    a = pmcmc_varied_parameter("a", c("p1", "p2"), 1:2),
    b = pmcmc_varied_parameter("b", c("p1", "p2"), 3:4),
    c = pmcmc_parameter("c", 5),
    d = pmcmc_parameter("d", 6))
  proposal_fixed <- diag(2)
  proposal_varied <- diag(2) + 1
  p <- pmcmc_parameters_nested$new(parameters, proposal_varied, proposal_fixed)

  expect_equal(p$model(p$initial()),
               unname(apply(p$initial(), 2, as.list)))

  p <- pmcmc_parameters_nested$new(parameters, proposal_varied, proposal_fixed,
                                   transform = function(x) as.list(log(x)))
  expect_equal(p$model(p$initial()),
               unname(apply(p$initial(), 2, function(x) as.list(log(x)))))
})


test_that("pmcmc_parameters_nested propose", {
  parameters <- list(
    a = pmcmc_varied_parameter("a", c("p1", "p2", "p3"), 1:3),
    b = pmcmc_parameter("b", 5, prior = function(x) 5),
    c = pmcmc_varied_parameter("c", c("p1", "p2", "p3"), 3:5),
    d = pmcmc_parameter("d", 6, prior = function(x) 6))
  proposal_fixed <- diag(2)
  proposal_varied <- diag(2) + 1
  p <- pmcmc_parameters_nested$new(parameters, proposal_varied, proposal_fixed)

  init <- p$initial()

  prop <- p$propose(init, type = "fixed")
  expect_identical(dimnames(prop), dimnames(init))
  expect_false(identical(prop[c(2, 4), ], init[c(2, 4), ]))
  expect_true(identical(prop[c(1, 3), ], init[c(1, 3), ]))
  expect_identical(prop[c(2, 4), 1], prop[c(2, 4), 2])

  prop <- p$propose(init, type = "varied")
  expect_identical(dimnames(prop), dimnames(init))
  expect_true(identical(prop[c(2, 4), ], init[c(2, 4), ]))
  expect_false(identical(prop[c(1, 3), ], init[c(1, 3), ]))

  prop <- p$propose(init, type = "both")
  expect_identical(dimnames(prop), dimnames(init))
  expect_false(any(identical(prop, init)))
  expect_identical(prop[c(2, 4), 1], prop[c(2, 4), 2])
})

test_that("pmcmc_parameters_nested propose - 1 varied 1 pop", {
  parameters <- list(a = pmcmc_varied_parameter("a", "p1", 1, prior = dnorm))
  proposal_varied <- diag(1)
  p <- pmcmc_parameters_nested$new(parameters, proposal_varied)
  init <- p$initial()

  prop <- p$propose(init, type = "varied")
  expect_true(inherits(prop, "matrix"))
  expect_identical(dimnames(prop), dimnames(init))

  expect_identical(p$propose(init, type = "fixed"), init)
  expect_equal(
    withr::with_seed(1, p$propose(init, type = "both")),
    withr::with_seed(1, p$propose(init, type = "varied")))
})


test_that("pmcmc_parameters_nested propose - fixed only", {
  parameters <- list(b = pmcmc_parameter("b", 1))
  proposal_varied <- NULL
  proposal_fixed <- diag(1)
  p <- pmcmc_parameters_nested$new(parameters, NULL, proposal_fixed,
                                   populations = letters[1:2])
  init <- p$initial()

  prop <- withr::with_seed(1, p$propose(init, type = "fixed"))
  expect_true(identical(prop[1, 1], prop[1, 2]))
  expect_true(all(prop != init))

  expect_equal(withr::with_seed(1, p$propose(init, type = "both")),
               prop)

  expect_identical(p$propose(init, type = "varied"), init)
})


test_that("pmcmc_parameters_nested fix errors", {
  parameters <- list(
    a = pmcmc_parameter("a", 1, prior = dexp),
    b = pmcmc_parameter("b", 1, prior = dexp))
  p <- pmcmc_parameters_nested$new(parameters, NULL, diag(2), c("x", "y"))
  expect_error(
    p$fix(matrix(1)),
    "colnames of 'fixed' must be identical to '$populations()'",
    fixed = TRUE)
  expect_error(
    p$fix(matrix(1, ncol = 2, dimnames = list(NULL, c("x", "y")))),
    "'fixed' must have rownames (parameters)",
    fixed = TRUE)
  expect_error(
    p$fix(matrix(1, ncol = 2, dimnames = list("x", c("x", "y")))),
    "Fixed parameters not found in model: 'x'")
  expect_error(
    p$fix(matrix(1, 2, 2, dimnames = list(c("b", "b"), c("x", "y")))),
    "Duplicate fixed parameters")
  expect_error(
    p$fix(matrix(1, 2, 2, dimnames = list(c("a", "b"), c("x", "y")))),
    "Cannot fix all parameters")
  expect_error(
    p$fix(matrix(1:2, 1, 2, dimnames = list("b", c("x", "y")))),
    "Fixed fixed parameters are not everywhere fixed")
})


test_that("can fix a subset of parameters", {
  parameters <- list(
    a = pmcmc_varied_parameter("a", c("p1", "p2"), 1:2,
                               prior = list(function(x) 1, function(x) 2)),
    b = pmcmc_parameter("b", 3,
                        prior = function(x) 1),
    c = pmcmc_varied_parameter("c", c("p1", "p2"), 4:5,
                               prior = list(function(x) 1, function(x) 2)),
    d = pmcmc_parameter("d", 5, prior = function(x) 6),
    e = pmcmc_varied_parameter("e", c("p1", "p2"), 7:8,
                               prior = list(function(x) 3, function(x) 4)),
    f = pmcmc_parameter("f", 6, prior = function(x) 9))
  proposal_fixed <- diag(3)
  proposal_varied <- diag(3) + 1
  p <- pmcmc_parameters_nested$new(parameters, proposal_varied, proposal_fixed)

  fixed <- cbind(p1 = c(a = 10, b = 13, c = 14),
                 p2 = c(a = 11, b = 13, c = 15))
  p2 <- p$fix(fixed)

  expect_equal(p2$names(), c("d", "e", "f"))
  expect_equal(p2$names("fixed"), c("d", "f"))
  expect_equal(p2$names("varied"), "e")
  expect_equal(p2$populations(), c("p1", "p2"))

  keep <- p2$names()
  expect_equal(p2$initial(), p$initial()[keep, ])

  cmp <- p$initial()
  cmp[rownames(fixed), ] <- fixed
  expect_identical(p2$model(p2$initial()), p$model(cmp))

  cmp <- subset(p$summary(), name %in% keep)
  rownames(cmp) <- NULL
  expect_equal(p2$summary(), cmp)

  init <- p2$initial()
  cmp <- pmcmc_parameters_nested$new(parameters[keep], diag(1) + 1, diag(2))
  expect_identical(
    withr::with_seed(1, p2$propose(init, "both")),
    withr::with_seed(1, cmp$propose(init, "both")))
  expect_identical(p2$prior(init), cmp$prior(init))
  expect_identical(p2$summary(), cmp$summary())

})


test_that("can fix all varied parameters", {
  parameters <- list(
    a = pmcmc_varied_parameter("a", "p1", 1),
    b = pmcmc_parameter("b", 2))
  proposal_fixed <- diag(1)
  proposal_varied <- diag(1) + 1
  p <- pmcmc_parameters_nested$new(parameters, proposal_varied, proposal_fixed)

  fix_p <- p$fix(matrix(1, dimnames = list("a", "p1")))

  cmp <- pmcmc_parameters_nested$new(parameters["b"], NULL, diag(1), "p1")

  expect_equal(
    withr::with_seed(1, fix_p$propose(fix_p$initial(), type = "both")),
    withr::with_seed(1, cmp$propose(cmp$initial(), type = "both")))

  expect_identical(fix_p$model(fix_p$initial()),
                   list(list(a = 1, b = 2)))
})


test_that("can fix all fixed parameters", {
  parameters <- list(
    a = pmcmc_varied_parameter("a", "p1", 1),
    b = pmcmc_parameter("b", 2))
  proposal_fixed <- diag(1)
  proposal_varied <- diag(1) + 1
  p <- pmcmc_parameters_nested$new(parameters, proposal_varied, proposal_fixed)

  fix_p <- p$fix(matrix(3, dimnames = list("b", "p1")))

  cmp <- pmcmc_parameters_nested$new(parameters["a"], diag(1) + 1, NULL)

  expect_equal(withr::with_seed(1, fix_p$propose(fix_p$initial(), "both")),
               withr::with_seed(1, cmp$propose(cmp$initial(), "both")))

  expect_identical(fix_p$model(fix_p$initial()),
                   list(list(a = 1, b = 3)))
})


test_that("Prevent use of variable fixed parameters", {
  pars <- list(
    pmcmc_varied_parameter("a", c("x", "y"), 1:2),
    pmcmc_varied_parameter("b", c("x", "y"), 3:4),
    pmcmc_parameter("c", 5),
    pmcmc_parameter("d", 6))
  p <- pmcmc_parameters_nested$new(pars, diag(2), diag(2))
  expect_error(
    p$validate(matrix(1:8, 4, 2)),
    "Fixed parameters are not everywhere fixed")
})


test_that("Control over transform", {
  pars <- list(
    pmcmc_varied_parameter("a", c("x", "y"), 1:2),
    pmcmc_varied_parameter("b", c("x", "y"), 3:4),
    pmcmc_parameter("c", 5),
    pmcmc_parameter("d", 6))
  p <- pmcmc_parameters_nested$new(pars, diag(2), diag(2), transform = as.list)
  expect_equal(r6_private(p)$transform, list(x = as.list, y = as.list))

  transform <- list(x = function(...) NULL, y = function(...) list())
  p <- pmcmc_parameters_nested$new(pars, diag(2), diag(2),
                                   transform = transform)
  expect_equal(r6_private(p)$transform, transform)

  expect_error(
    pmcmc_parameters_nested$new(pars, diag(2), diag(2),
                                transform = unname(transform)),
    "If 'transform' is a list, its names must be the populations")
})
