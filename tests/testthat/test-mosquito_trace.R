test_that("null mosquito trace works", {
  p <- 1
  tmax <- 10
  mod <- make_MicroMoB(tmax = tmax, p = p)
  expect_error(setup_mosquito_trace(model = mod, oviposit = rpois(10, 10)))
  expect_error(setup_mosquito_trace(model = mod, oviposit = NaN))
  expect_error(setup_mosquito_trace(model = mod, oviposit = Inf))

  setup_mosquito_trace(model = mod, oviposit = 5)
  expect_equal(compute_oviposit(mod), 5)

  p <- 3
  tmax <- 10
  mod <- make_MicroMoB(tmax = tmax, p = p)
  expect_error(setup_mosquito_trace(model = mod, oviposit = 5))
  setup_mosquito_trace(model = mod, oviposit = c(1,5,10))
  expect_equal(compute_oviposit(mod), c(1,5,10))
})
