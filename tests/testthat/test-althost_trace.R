test_that("trace alternative blood host model works", {

  tmax <- 20
  p <- 4

  mod <- make_MicroMoB(tmax = tmax, p = p)
  setup_alternative_trace(mod)

  expect_equal(mod$alternative$O, matrix(0, nrow = p, ncol = tmax))

  expect_equal(compute_O(mod), rep(0, p))

  mod <- make_MicroMoB(tmax = tmax, p = p)
  setup_alternative_trace(mod, O = c(1, 2, 3, 4))

  expect_equal(mod$alternative$O, replicate(tmax, c(1, 2, 3, 4)))

  expect_equal(compute_O(mod), c(1, 2, 3, 4))

  mod <- make_MicroMoB(tmax = tmax, p = p)
  setup_alternative_trace(mod, O = matrix(1:(p*tmax), nrow = p, ncol = tmax))

  expect_equal(compute_O(mod), c(1, 2, 3, 4))

})
