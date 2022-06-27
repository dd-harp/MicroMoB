test_that("human object setup is working", {

  n <- 3
  p <- 2
  tmax <- 10

  mod <- make_MicroMoB(tmax = tmax, p = p)

  theta <- matrix(rexp(n*p), n, p)
  theta <- theta / rowSums(theta)
  wf <- rep(1, n)
  SIP <- matrix(data = rpois(n = n*3, lambda = 50), nrow = n, ncol = 3)
  c <- 0.15

  setup_humans_SIP(model = mod, stochastic = FALSE, theta = theta, SIP = SIP)
  expect_equal(compute_x(mod), (SIP[, 2L] / rowSums(SIP)) * c)
  expect_equal(compute_wf(mod), wf)
  expect_equal(compute_Psi(mod), theta)
  expect_equal(compute_H(mod), rowSums(SIP))

})


test_that("deterministic updates of human SIP model work", {

  n <- 3
  p <- 2
  tmax <- 100

  mod <- make_MicroMoB(tmax = tmax, p = p)

  theta <- matrix(rexp(n*p), n, p)
  theta <- theta / rowSums(theta)
  wf <- rep(1, n)
  SIP <- matrix(data = rpois(n = n*3, lambda = 50), nrow = n, ncol = 3)
  c <- 0.15

  setup_humans_SIP(model = mod, stochastic = FALSE, theta = theta, SIP = SIP)

  for (i in 1:tmax) {
    step_humans(model = mod)
    mod$global$tnow <- mod$global$tnow + 1L
  }

  expect_true(all(mod$human$SIP[, "S"] > SIP[, 1L]))
  expect_true(all(mod$human$SIP[, "I"] < SIP[, 2L]))
  expect_true(all(mod$human$SIP[, "P"] < SIP[, 3L]))

})

test_that("stochastic updates of human SIP model work", {

  n <- 3
  p <- 2
  tmax <- 100

  mod <- make_MicroMoB(tmax = tmax, p = p)

  theta <- matrix(rexp(n*p), n, p)
  theta <- theta / rowSums(theta)
  wf <- rep(1, n)
  SIP <- matrix(data = rpois(n = n*3, lambda = 50), nrow = n, ncol = 3)
  c <- 0.15

  setup_humans_SIP(model = mod, stochastic = TRUE, theta = theta, SIP = SIP)

  for (i in 1:tmax) {
    step_humans(model = mod)
    mod$global$tnow <- mod$global$tnow + 1L
  }

  expect_true(all(mod$human$SIP[, "S"] >= SIP[, 1L]))
  expect_true(all(mod$human$SIP[, "I"] <= SIP[, 2L]))
  expect_true(all(mod$human$SIP[, "P"] <= SIP[, 3L]))

})
