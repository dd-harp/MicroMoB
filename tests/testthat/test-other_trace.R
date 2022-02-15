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


test_that("test JSON config working", {

  library(jsonlite)

  # sending to JSON does not change R type when read back in
  par <- list(
    "O" = rep(1, 5)
  )


  json_path <- tempfile(pattern = "visitor_par", fileext = ".json")
  write_json(x = par, path = json_path, digits = NA)
  par_in <- get_config_alternative_trace(path = json_path)
  expect_true(all.equal(par, par_in))

  # reject obviously bad input
  par <- list(
    "O" = NULL
  )

  json_path <- tempfile(pattern = "visitor_par", fileext = ".json")
  write_json(x = par, path = json_path, digits = NA)
  expect_error(get_config_alternative_trace(path = json_path))

  unlink(x = json_path)

})


test_that("JSON parameters can read in", {
  path <- system.file("extdata", "other_trace.json", package = "MicroMoB")
  pars <- get_config_alternative_trace(path = path)

  expect_true(length(pars) == 1L)

  expect_true(is.numeric(pars$O))
  expect_true(is.vector(pars$O) | is.matrix(pars$O))
})
