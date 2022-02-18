test_that("null mosquito trace works", {
  p <- 1
  tmax <- 10
  mod <- make_MicroMoB(tmax = tmax, p = p)
  expect_error(setup_mosquito_trace(model = mod, oviposit = rpois(10, 10)))
  expect_error(setup_mosquito_trace(model = mod, oviposit = NaN))
  expect_error(setup_mosquito_trace(model = mod, oviposit = Inf))

  setup_mosquito_trace(model = mod, oviposit = 5)
  expect_equal(compute_oviposit(mod), 5)

  expect_error(step_mosquitoes(mod))
  expect_error(compute_f(mod))
  expect_error(compute_q(mod))
  expect_error(compute_Z(mod))

  expect_equal(output_mosquitoes(mod), data.frame())

  p <- 3
  tmax <- 10
  mod <- make_MicroMoB(tmax = tmax, p = p)
  expect_error(setup_mosquito_trace(model = mod, oviposit = 5))
  setup_mosquito_trace(model = mod, oviposit = c(1,5,10))
  expect_equal(compute_oviposit(mod), c(1,5,10))

  expect_error(step_mosquitoes(mod))
  expect_error(compute_f(mod))
  expect_error(compute_q(mod))
  expect_error(compute_Z(mod))

  expect_equal(output_mosquitoes(mod), data.frame())

})


test_that("test JSON config working", {

  library(jsonlite)

  # sending to JSON does not change R type when read back in
  par <- list(
    "oviposit" = rep(1, 3)
  )

  json_path <- tempfile(pattern = "mosquito_par", fileext = ".json")
  write_json(x = par, path = json_path, digits = NA)
  par_in <- get_config_mosquito_trace(path = json_path)
  expect_true(all.equal(par, par_in))

  # reject obviously bad input
  par <- list(
    "oviposit" = NULL
  )

  json_path <- tempfile(pattern = "mosquito_par", fileext = ".json")
  write_json(x = par, path = json_path, digits = NA)
  expect_error(get_config_mosquito_trace(path = json_path))

  unlink(x = json_path)

})


test_that("JSON parameters can read in", {
  path <- system.file("extdata", "mosquito_trace.json", package = "MicroMoB")
  pars <- get_config_mosquito_trace(path = path)

  expect_true(length(pars) == 1L)

  expect_true(is.numeric(pars$oviposit))
  expect_true(is.vector(pars$oviposit))

})
