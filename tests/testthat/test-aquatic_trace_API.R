test_that("test aqua trace parameter API", {
  p <- 3
  tmax <- 20
  lambda <- 50
  mod <- make_MicroMoB(tmax = tmax, p = p)
  setup_aqua_trace(model = mod, lambda = lambda, stochastic = FALSE)

  # change lambda
  set_lambda_aqua_trace(model = mod, lambda = 0.2, times = 2, places = 1)
  expect_equal(get_lambda_aqua_trace(model = mod, times = 2, places = 1), 0.2)
  expect_true(all(as.vector(mod$aqua$lambda)[-4] == lambda))

  set_lambda_aqua_trace(model = mod, lambda = 0.4, times = 5:6, places = 2:3)
  expect_true(all(get_lambda_aqua_trace(model = mod, times = 5:6, places = 2:3) == 0.4))
  expect_equal(sum(mod$aqua$lambda == lambda), 55)

  set_lambda_aqua_trace(model = mod, lambda = matrix(c(0.1,0.2,0.3,0.4),2,2), times = c(10, 20), places = c(1, 3))
  expect_true(all(get_lambda_aqua_trace(model = mod, times = c(10,20),places = c(1,3)) == c(0.1,0.2,0.3,0.4)))
  expect_equal(sum(mod$aqua$lambda == lambda), 51)

  expect_error(set_lambda_aqua_trace(model = mod, lambda = -0.5, times = 1, places = 1))
  expect_error(set_lambda_aqua_trace(model = mod, lambda = 0.5, times = -5, places = 1))
  expect_error(set_lambda_aqua_trace(model = mod, lambda = 0.5, times = 5, places = -1))
  expect_error(set_lambda_aqua_trace(model = mod, lambda = 0.5, times = 5, places = 100))
  expect_error(set_lambda_aqua_trace(model = mod, lambda = 0.5, times = 100, places = 1))

})
