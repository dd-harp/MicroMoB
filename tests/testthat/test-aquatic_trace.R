test_that("test emergence model of aquatic dynamics errors with incorrect lambda", {
  p <- 2
  tmax <- 10
  mod <- make_MicroMoB(tmax = tmax, p = p)

  lambda <- c(10, 100, 1000)
  expect_error(setup_aqua_trace(model = mod, lambda = lambda, stochastic = FALSE))
  expect_error(setup_aqua_trace(model = mod, lambda = lambda, stochastic = TRUE))

  lambda <- matrix()
  expect_error(setup_aqua_trace(model = mod, lambda = lambda, stochastic = FALSE))
  expect_error(setup_aqua_trace(model = mod, lambda = lambda, stochastic = TRUE))

  lambda <- matrix(rexp(100), 10, 10)
  expect_error(setup_aqua_trace(model = mod, lambda = lambda, stochastic = FALSE))
  expect_error(setup_aqua_trace(model = mod, lambda = lambda, stochastic = TRUE))
})


test_that("test emergence model of aquatic dynamics with vector lambda", {
  tmax <- 10
  p <- 2
  lambda <- c(10, 100)

  # deterministic
  mod <- make_MicroMoB(tmax = tmax, p = p)
  setup_aqua_trace(model = mod, lambda = lambda, stochastic = FALSE)
  expect_equal(compute_emergents(model = mod), lambda)

  # stochastic
  mod <- make_MicroMoB(tmax = tmax, p = p)
  setup_aqua_trace(model = mod, lambda = lambda, stochastic = TRUE)
  lambda <- compute_emergents(model = mod)
  expect_true(lambda[1] < lambda[2])

})


test_that("test emergence model of aquatic dynamics with 365 matrix lambda, tmax < 365", {
  tmax <- 50
  p <- 2
  lambda <- matrix(rnorm(n = 365 * p, mean = rep(c(10, 100), each = 365), sd = rep(c(2.5, 20), each = 365)), nrow = p, ncol = 365, byrow = TRUE)
  lambda <- pmax(lambda, 0)

  # deterministic
  mod <- make_MicroMoB(tmax = tmax, p = p)
  setup_aqua_trace(model = mod, lambda = lambda, stochastic = FALSE)

  expect_equal(mod$aqua$lambda, lambda[, 1:tmax])
  expect_equal(compute_emergents(model = mod), lambda[, 1])

  # stochastic
  mod <- make_MicroMoB(tmax = tmax, p = p)
  setup_aqua_trace(model = mod, lambda = lambda, stochastic = TRUE)

  expect_equal(mod$aqua$lambda, lambda[, 1:tmax])
  lambda <- compute_emergents(model = mod)
  expect_true(lambda[1] < lambda[2])

})


test_that("test emergence model of aquatic dynamics with 365 matrix lambda, tmax > 365", {
  tmax <- 730
  p <- 2
  lambda <- matrix(rnorm(n = 365 * p, mean = rep(c(10, 100), each = 365), sd = rep(c(2.5, 20), each = 365)), nrow = p, ncol = 365, byrow = TRUE)
  lambda <- pmax(lambda, 0)

  # deterministic
  mod <- make_MicroMoB(tmax = tmax, p = p)
  setup_aqua_trace(model = mod, lambda = lambda, stochastic = FALSE)

  expect_equal(mod$aqua$lambda, cbind(lambda, lambda))
  expect_equal(compute_emergents(model = mod), lambda[, 1])

  mod$global$tnow <- 366
  expect_equal(compute_emergents(model = mod), lambda[, 1])

  mod$global$tnow <- tmax
  expect_equal(compute_emergents(model = mod), lambda[, 365])

  # stochastic
  mod <- make_MicroMoB(tmax = tmax, p = p)
  setup_aqua_trace(model = mod, lambda = lambda, stochastic = TRUE)

  expect_equal(mod$aqua$lambda, cbind(lambda, lambda))
  lambda <- compute_emergents(model = mod)
  expect_true(lambda[1] < lambda[2])

})


test_that("test emergence model of aquatic dynamics with tmax matrix lambda", {
  tmax <- 20
  p <- 2
  lambda <- matrix(rnorm(n = 20 * p, mean = rep(c(10, 100), each = 20), sd = rep(c(2.5, 20), each = 20)), nrow = p, ncol = 20, byrow = TRUE)
  lambda <- pmax(lambda, 0)

  # deterministic
  mod <- make_MicroMoB(tmax = tmax, p = p)
  setup_aqua_trace(model = mod, lambda = lambda, stochastic = FALSE)

  expect_equal(mod$aqua$lambda, lambda)
  expect_equal(compute_emergents(model = mod), lambda[, 1])

  mod$global$tnow <- 19
  expect_equal(compute_emergents(model = mod), lambda[, 19])

  # stochastic
  mod <- make_MicroMoB(tmax = tmax, p = p)
  setup_aqua_trace(model = mod, lambda = lambda, stochastic = TRUE)

  expect_equal(mod$aqua$lambda, lambda)
  lambda <- compute_emergents(model = mod)
  expect_true(lambda[1] < lambda[2])

})


test_that("test JSON config working", {

  library(jsonlite)

  t <- 50
  p <- 3

  # scalar lambda
  par <- list(
    "stochastic" = FALSE,
    "lambda" = 5
  )

  json_path <- tempfile(pattern = "aqua_par", fileext = ".json")
  write_json(x = par, path = json_path)
  par_in <- get_config_aqua_trace(path = json_path)
  expect_true(all.equal(par, par_in))

  # vector lambda
  par <- list(
    "stochastic" = FALSE,
    "lambda" = rpois(n = t, lambda = 10)
  )

  json_path <- tempfile(pattern = "aqua_par", fileext = ".json")
  write_json(x = par, path = json_path)
  par_in <- get_config_aqua_trace(path = json_path)
  expect_true(all.equal(par, par_in))

  # matrix lambda
  par <- list(
    "stochastic" = FALSE,
    "lambda" = matrix(rpois(n = t * p, lambda = 10), nrow = p, ncol = t)
  )

  json_path <- tempfile(pattern = "aqua_par", fileext = ".json")
  write_json(x = par, path = json_path)
  par_in <- get_config_aqua_trace(path = json_path)
  expect_true(all.equal(par, par_in))

  # reject obviously bad input
  par <- list(
    "stochastic" = FALSE,
    "lambda" = array(data = 1:27, dim = c(3, 3, 3))
  )

  json_path <- tempfile(pattern = "aqua_par", fileext = ".json")
  write_json(x = par, path = json_path)
  expect_error(get_config_aqua_trace(path = json_path))

  par <- list(
    "stochastic" = FALSE,
    "lambda" = "c(1, 2, 3)"
  )

  json_path <- tempfile(pattern = "aqua_par", fileext = ".json")
  write_json(x = par, path = json_path)
  expect_error(get_config_aqua_trace(path = json_path))

  par <- list(
    "stochastic" = 5,
    "lambda" = c(1, 2, 3)
  )

  json_path <- tempfile(pattern = "aqua_par", fileext = ".json")
  write_json(x = par, path = json_path)
  expect_error(get_config_aqua_trace(path = json_path))

  unlink(x = json_path)

})

test_that("JSON parameters can read in", {
  path <- system.file("extdata", "aqua_trace.json", package = "MicroMoB")
  pars <- get_config_aqua_trace(path = path)
  expect_true(length(pars) == 2L)

  expect_true(is.logical(pars$stochastic))
  expect_true(length(pars$stochastic) == 1L)

  expect_true(is.numeric(pars$lambda))
  expect_true(is.vector(pars$lambda) | is.matrix(pars$lambda))

})
