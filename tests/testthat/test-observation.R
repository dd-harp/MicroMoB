test_that("pfpr observation for SIP works", {
  n <- 3
  p <- 2
  tmax <- 100

  mod <- make_MicroMoB(tmax = tmax, p = p)

  theta <- matrix(1, n, p)
  theta <- theta / rowSums(theta)
  SIP <- matrix(data = rpois(n = n*3, lambda = 50), nrow = n, ncol = 3)

  setup_humans_SIP(model = mod, stochastic = FALSE, theta = theta, SIP = SIP)

  # 100% test cov, sens/spec 1
  test_pars <- list(testprop = c(1, 1, 1), sens = 1, spec = 1)
  pfpr <- observe_pfpr(model = mod, parameters = test_pars)

  expect_true(all(pfpr["neg", "pos", ] == 0))
  expect_true(all(pfpr["pos", "neg", ] == 0))
  expect_true(all(pfpr["pos", "pos", ] == SIP[, 2L]))
  expect_true(all(pfpr["neg", "neg", ] == rowSums(SIP[, c(1, 3)])))

  # 100%t test cov, sens 0, spec 0
  test_pars <- list(testprop = c(1, 1, 1), sens = 0, spec = 0)
  pfpr <- observe_pfpr(model = mod, parameters = test_pars)

  expect_true(all(pfpr["pos", "pos", ] == 0))
  expect_true(all(pfpr["neg", "neg", ] == 0))
  expect_true(all(pfpr["pos", "neg", ] == SIP[, 2L]))
  expect_true(all(pfpr["neg", "pos", ] == rowSums(SIP[, c(1, 3)])))

  # 0% test coverage, nobody should be tested
  test_pars <- list(testprop = c(0, 0, 0), sens = 0, spec = 0)
  pfpr <- observe_pfpr(model = mod, parameters = test_pars)
  expect_true(sum(pfpr) == 0)

})


test_that("pfpr observation for SIS works", {
  n <- 3
  p <- 2
  tmax <- 100

  mod <- make_MicroMoB(tmax = tmax, p = p)

  theta <- matrix(1, n, p)
  theta <- theta / rowSums(theta)
  H <- rpois(n = n, lambda = 100)
  X <- rpois(n = n, lambda = 15)

  setup_humans_SIS(model = mod, stochastic = FALSE, theta = theta, H = H, X = X)

  # 100% test cov, sens/spec 1
  test_pars <- list(testprop = c(1, 1, 1), sens = 1, spec = 1)
  pfpr <- observe_pfpr(model = mod, parameters = test_pars)

  expect_true(all(pfpr["neg", "pos", ] == 0))
  expect_true(all(pfpr["pos", "neg", ] == 0))
  expect_true(all(pfpr["pos", "pos", ] == X))
  expect_true(all(pfpr["neg", "neg", ] == H - X))

  # 100%t test cov, sens 0, spec 0
  test_pars <- list(testprop = c(1, 1, 1), sens = 0, spec = 0)
  pfpr <- observe_pfpr(model = mod, parameters = test_pars)

  expect_true(all(pfpr["pos", "pos", ] == 0))
  expect_true(all(pfpr["neg", "neg", ] == 0))
  expect_true(all(pfpr["pos", "neg", ] == X))
  expect_true(all(pfpr["neg", "pos", ] == H - X))

  # 0% test coverage, nobody should be tested
  test_pars <- list(testprop = c(0, 0, 0), sens = 0, spec = 0)
  pfpr <- observe_pfpr(model = mod, parameters = test_pars)
  expect_true(sum(pfpr) == 0)

})
