test_that("constructing human objects (strata) works when J is not specified", {

  expect_error(setup.human("strata", model = new.env(), H = NULL))
  expect_error(setup.human("strata", model = new.env(), H = rep(Inf, 5)))
  expect_error(setup.human("strata", model = new.env(), H = rep(-5, 5)))
  expect_error(setup.human("strata", model = new.env(), H = rep(NaN, 5)))
  expect_error(setup.human("strata", model = new.env(), H = rep(NA, 5)))

  expect_error(setup.human("strata", model = new.env(), H = c(Inf, 5)))
  expect_error(setup.human("strata", model = new.env(), H = c(-5, 5)))
  expect_error(setup.human("strata", model = new.env(), H = c(NaN, 5)))
  expect_error(setup.human("strata", model = new.env(), H = c(NA, 5)))

  expect_error(setup.human("strata", model = new.env(), H = c(1, 5), J = diag(5)))

  model <- new.env()
  human <- setup.human("strata", model = model, H = c(1, 5))
  expect_equal(model$human$J, diag(2))
  expect_equal(model$human$H, c(1, 5))
  expect_true(inherits(model$human, "strata"))
})


test_that("constructing human objects (strata) works when specifying J", {

  J <- matrix(
    c(0.3, 0.5, 0.2,
      0.1, 0.6, 0.3), nrow = 3, ncol = 2, byrow = F
  )
  H <- c(s1 = 50, s2 = 60)

  bad_values <- c(Inf, NaN, NA, -5)
  for (v in bad_values) {
    J[1,2] <- v
    expect_error(setup.human("strata", model = new.env(), H = H, J = J))
    expect_error(setup.human("strata", model = new.env(), H = H, J = J))
    expect_error(setup.human("strata", model = new.env(), H = H, J = J))
    expect_error(setup.human("strata", model = new.env(), H = H, J = J))
  }

  J[1,2] <- 0.1

  expect_error(setup.human("simple", model = new.env(), H = H, J = J))

  expect_error(setup.human("strata", model = new.env(), H = NULL, J = J))
  expect_error(setup.human("strata", model = new.env(), H = rep(Inf, 5), J = J))
  expect_error(setup.human("strata", model = new.env(), H = rep(-5, 5), J = J))
  expect_error(setup.human("strata", model = new.env(), H = rep(NaN, 5), J = J))
  expect_error(setup.human("strata", model = new.env(), H = rep(NA, 5), J = J))

  expect_error(setup.human("strata", model = new.env(), H = c(Inf, 5), J = J))
  expect_error(setup.human("strata", model = new.env(), H = c(-5, 5), J = J))
  expect_error(setup.human("strata", model = new.env(), H = c(NaN, 5), J = J))
  expect_error(setup.human("strata", model = new.env(), H = c(NA, 5), J = J))

  model <- new.env()
  human <- setup.human("strata", model = model, H = H, J = J)
  expect_true(sum(model$human$J %*% model$human$H) == sum(H))

  human_pop <- matrix(
    c(15, 6,
      25, 36,
      10, 18), nrow = 3, ncol = 2, byrow = TRUE
  )
  expect_equal(as.vector(model$human$J %*% model$human$H), rowSums(human_pop))
})
