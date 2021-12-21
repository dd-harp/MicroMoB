test_that("trace visitor model works", {

  tmax <- 20
  p <- 4

  mod <- make_MicroMoB(tmax = tmax, p = p)
  setup_visitor_trace(mod)

  expect_equal(mod$visitor$Wd, matrix(0, nrow = p, ncol = tmax))
  expect_equal(mod$visitor$xd, matrix(0, nrow = p, ncol = tmax))

  expect_equal(compute_Wd(mod), rep(0, p))
  expect_equal(compute_xd(mod), rep(0, p))

  mod <- make_MicroMoB(tmax = tmax, p = p)
  setup_visitor_trace(mod, Wd = c(1, 2, 3, 4), xd = c(0.1, 0.2, 0.3, 0.4))

  expect_equal(mod$visitor$Wd, replicate(tmax, c(1, 2, 3, 4)))
  expect_equal(mod$visitor$xd, replicate(tmax, c(0.1, 0.2, 0.3, 0.4)))

  expect_equal(compute_Wd(mod), c(1, 2, 3, 4))
  expect_equal(compute_xd(mod), c(0.1, 0.2, 0.3, 0.4))

  mod <- make_MicroMoB(tmax = tmax, p = p)
  setup_visitor_trace(mod, Wd = matrix(1:(p*tmax), nrow = p, ncol = tmax), xd = matrix(1:(p*tmax) / (p*tmax), nrow = p, ncol = tmax))

  expect_equal(compute_Wd(mod), c(1, 2, 3, 4))
  expect_equal(compute_xd(mod), c(0.0125, 0.0250, 0.0375, 0.0500))

})

