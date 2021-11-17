test_that("setting up human objects (strata) works when J is not specified", {

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
  setup.human("strata", model = model, H = c(1, 5))
  expect_equal(model$human$J, diag(2))
  expect_equal(model$human$H, c(1, 5))
  expect_true(inherits(model$human, "strata"))
})


test_that("setting up human objects (strata) works when specifying J", {

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
  setup.human("strata", model = model, H = H, J = J)
  expect_true(sum(model$human$J %*% model$human$H) == sum(H))

  human_pop <- matrix(
    c(15, 6,
      25, 36,
      10, 18), nrow = 3, ncol = 2, byrow = TRUE
  )
  expect_equal(as.vector(model$human$J %*% model$human$H), rowSums(human_pop))

  # test it works when we give J already formulated as binary
  model1 <- new.env()
  setup.human("strata", model = model1, H = model$human$H, J = model$human$J)
  expect_true(all.equal(as.list(model), as.list(model1)))
})


test_that("setting up time spent (day) works", {

  model <- new.env()
  setup.human("strata", model = model, H = rep(10, 3))
  setup.timespent("day", model = model)

  expect_true(inherits(model$tisp, "day"))
  expect_equal(model$tisp$theta, diag(3))
})

test_that("setting up time spent (dt) works", {

  # 2 patch, 3 strata
  #   s1p1  s1p2 s2p1 s2p2
  H <- c(50, 20, 10, 30, 20, 10)
  J <- matrix(
    c(1, 0, 1, 0, 1, 0,
      0, 1, 0, 1, 0, 1),
    nrow = 2, ncol = 6, byrow = TRUE
  )

  expect_equal(as.vector(J %*% H), c(sum(H[1], H[3], H[5]), sum(H[2], H[4], H[6])))

  theta_day <- matrix(
    c(0.4, 0.6,
      0.5, 0.5,
      0.2, 0.8,
      0.5, 0.5,
      0.6, 0.4,
      0.2, 0.8), nrow = length(H), ncol = 2, byrow = TRUE
  )

  theta_night <- matrix(
    c(0.9, 0.1,
      0.6, 0.4,
      0.5, 0.5,
      0.6, 0.4,
      0.8, 0.2,
      0.4, 0.6), nrow = length(H), ncol = 2, byrow = TRUE
  )

  theta <- list(theta_day, theta_night)

  model <- new.env()
  setup.human("strata", model = model, H = H, J = J)
  setup.timespent("dt", model = model, theta = theta)
  setup.biteweight("null", model = model)

  xi <- c(0.15, 0.85)
  W <- compute.timespent(tisp = model$tisp, biteweight = model$biteweight, human = model$human, xi = xi, t = 1)

  wt <- rep(1, length(H))
  W_expected <- ((t(theta_day) %*% (wt * H)) * (xi[1])) + ((t(theta_night) %*% (wt * H)) * (xi[2]))

  expect_equal(W, W_expected)
})
