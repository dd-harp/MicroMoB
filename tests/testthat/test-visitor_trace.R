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


test_that("test JSON config working", {

  library(jsonlite)

  # sending to JSON does not change R type when read back in
  par <- list(
    "Wd" = rep(1, 5),
    "xd" = rep(0.01, 365)
  )

  json_path <- tempfile(pattern = "visitor_par", fileext = ".json")
  write_json(x = par, path = json_path, digits = NA)
  par_in <- get_config_visitor_trace(path = json_path)
  expect_true(all.equal(par, par_in))

  # reject obviously bad input
  par <- list(
    "Wd" = rep(1, 5),
    "xd" = NULL
  )

  json_path <- tempfile(pattern = "visitor_par", fileext = ".json")
  write_json(x = par, path = json_path, digits = NA)
  expect_error(get_config_visitor_trace(path = json_path))

  unlink(x = json_path)

})
