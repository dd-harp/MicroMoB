test_that("human object setup is working", {

  n <- 3
  p <- 3
  tmax <- 10

  mod <- make_MicroMoB(tmax = tmax, p = p)

  b <- 0.55
  c <- 0.15
  gamma <- 1/5

  theta <- matrix(rexp(9), n, p)
  theta <- theta / rowSums(theta)
  wf <- rep(1, n)
  H <- c(100, 80, 50)
  SIR <- matrix(
    c(85, 5, 10,
      70, 5, 5,
      25, 10, 15),
    nrow = n, ncol = 3, byrow = TRUE
  )

  setup_humans_SIR(model = mod, stochastic = FALSE, theta = theta, wf = wf, H = H, SIR = SIR, b = b, c = c, gamma = gamma)
  expect_equal(compute_W(mod), as.vector(t(theta) %*% (wf * H)))
  expect_equal(compute_wf(mod), rep(1, p))
  expect_equal(compute_x(mod), (SIR[, 2] / H) * c)
  expect_equal(compute_H(mod), H)
  expect_equal(compute_Psi(mod), theta)

})


test_that("deterministic updates of human SIR model work", {

  tmax <- 1e2
  n <- 3
  p <- 3

  mod <- make_MicroMoB(tmax = tmax, p = p)

  b <- 0.55
  c <- 0.15
  gamma <- 1/5

  theta <- matrix(rexp(9), n, p)
  theta <- theta / rowSums(theta)
  wf <- rep(1, 3)
  H <- c(100, 80, 50)
  SIR <- matrix(
    c(85, 5, 10,
      70, 5, 5,
      25, 10, 15),
    nrow = n, ncol = 3, byrow = TRUE
  )

  setup_humans_SIR(model = mod, stochastic = FALSE, theta = theta, wf = wf, H = H, SIR = SIR, b = b, c = c, gamma = gamma)

  for (i in 1:tmax) {
    step_humans(model = mod)
    mod$global$tnow <- mod$global$tnow + 1L
  }

  expect_equal(mod$human$SIR[, 3], rowSums(SIR[, 2:3]))
  expect_equal(mod$human$SIR[, 1], SIR[, 1])
  expect_true(all(mod$human$SIR[, 2] >= 0))
})


test_that("deterministic updates of human SIR model work with pulsed h", {

  tmax <- 1e2
  n <- 3
  p <- 3

  mod <- make_MicroMoB(tmax = tmax, p = p)

  b <- 0.55
  c <- 0.15
  gamma <- 1/5

  theta <- matrix(rexp(9), n, p)
  theta <- theta / rowSums(theta)
  wf <- rep(1, 3)
  H <- c(100, 80, 50)
  SIR <- matrix(
    c(85, 5, 10,
      70, 5, 5,
      25, 10, 15),
    nrow = n, ncol = 3, byrow = TRUE
  )

  setup_humans_SIR(model = mod, stochastic = FALSE, theta = theta, wf = wf, H = H, SIR = SIR, b = b, c = c, gamma = gamma)

  h <- rep(qexp(p = 0.25), n)
  mod$human$EIR <- h / mod$human$b

  step_humans(model = mod)

  expect_true(all(mod$human$SIR[, 1] < SIR[, 1]))
  expect_true(all(mod$human$SIR[, 2] > SIR[, 2]))
  expect_true(all(mod$human$SIR[, 3] > SIR[, 3]))

  mod$global$tnow <- mod$global$tnow + 1L
  prev_I <- mod$human$SIR[, "I"]

  mod$human$EIR <- rep(0, n)
  for (i in 2:tmax) {
    step_humans(model = mod)
    expect_true(all(mod$human$SIR[, "I"] < prev_I))
    mod$global$tnow <- mod$global$tnow + 1L
    prev_I <- mod$human$SIR[, "I"]
  }

  expect_equal(mod$human$SIR[, 1], SIR[, 1] * (1 - 0.25))
  expect_true(all(mod$human$SIR[, 2] >= 0))
  expect_true(all(mod$human$SIR[, 3] >= 0))
  expect_true(all(mod$human$SIR[, 3] > mod$human$SIR[, 2]))
  expect_equal(sum(SIR), sum(mod$human$SIR))

})


test_that("stochastic updates of human SIR model work", {

  tmax <- 1e2
  n <- 3
  p <- 3

  mod <- make_MicroMoB(tmax = tmax, p = p)

  b <- 0.55
  c <- 0.15
  gamma <- 1/5

  theta <- matrix(rexp(9), n, p)
  theta <- theta / rowSums(theta)
  wf <- rep(1, 3)
  H <- c(100, 80, 50)
  SIR <- matrix(
    c(85, 5, 10,
      70, 5, 5,
      25, 10, 15),
    nrow = n, ncol = 3, byrow = TRUE
  )

  setup_humans_SIR(model = mod, stochastic = TRUE, theta = theta, wf = wf, H = H, SIR = SIR, b = b, c = c, gamma = gamma)

  for (i in 1:tmax) {
    step_humans(model = mod)
    mod$global$tnow <- mod$global$tnow + 1L
  }

  expect_equal(mod$human$SIR[, 1], SIR[, 1])
  expect_true(all(mod$human$SIR[, 2] == 0))
  expect_equal(mod$human$SIR[, 3], rowSums(SIR[, 2:3]))
})


test_that("stochastic updates of human SIR model work with pulsed h", {

  tmax <- 1e2
  n <- 3
  p <- 3

  mod <- make_MicroMoB(tmax = tmax, p = p)

  b <- 0.55
  c <- 0.15
  gamma <- 1/5

  theta <- matrix(rexp(9), n, p)
  theta <- theta / rowSums(theta)
  wf <- rep(1, 3)
  H <- c(100, 80, 50)
  SIR <- matrix(
    c(85, 5, 10,
      70, 5, 5,
      25, 10, 15),
    nrow = n, ncol = 3, byrow = TRUE
  )

  setup_humans_SIR(model = mod, stochastic = TRUE, theta = theta, wf = wf, H = H, SIR = SIR, b = b, c = c, gamma = gamma)

  h <- rep(qexp(p = 0.5), n)
  mod$human$EIR <- h / mod$human$b

  step_humans(model = mod)

  expect_true(all(mod$human$SIR[, 1] <= SIR[, 1]))
  expect_true(all(mod$human$SIR[, 2] >= SIR[, 2]))
  expect_true(all(mod$human$SIR[, 3] >= SIR[, 3]))

  mod$global$tnow <- mod$global$tnow + 1L
  prev_I <- mod$human$SIR[, "I"]

  mod$human$EIR <- rep(0, n)
  for (i in 2:tmax) {
    step_humans(model = mod)
    expect_true(all(mod$human$SIR[, "I"] <= prev_I))
    mod$global$tnow <- mod$global$tnow + 1L
    prev_I <- mod$human$SIR[, "I"]
  }

  expect_true(all(mod$human$SIR[, 2] >= 0))
  expect_true(all(mod$human$SIR[, 3] >= 0))
  expect_true(all(mod$human$SIR[, 3] > mod$human$SIR[, 2]))
  expect_equal(sum(SIR), sum(mod$human$SIR))

})
