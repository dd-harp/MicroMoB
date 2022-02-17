library(jsonlite)
library(plumber)
library(callr)
library(httr)

test_that("main JSON config loads", {

  par <- list(
    p = 5,
    tmax = 10,
    aqua_path = system.file("extdata", "aqua_BH.json", package = "MicroMoB"),
    aqua_model = "BH",
    adult_path = system.file("extdata", "mosquito_RM.json", package = "MicroMoB"),
    adult_model = "RM"
  )

  json_path <- tempfile(pattern = "config", fileext = ".json")
  jsonlite::write_json(x = par, path = json_path, digits = NA, pretty = TRUE)

  par_in <- get_config_mosquito_MicroMoB(path = json_path)
  expect_true(all.equal(par, par_in))

  unlink(x = json_path)

})


test_that("API can be started and runs a mosquito-only simulation", {

  patches <- 1
  tmax <- 20

  M <- 120
  p <- 0.9
  lambda <- M*(1-p)

  nu <- 25
  f <- 0.3
  q <- 0.9
  eggs <- nu * f * M

  # static pars
  molt <-  0.1
  surv <- 0.9

  # solve L
  L <- lambda * ((1/molt) - 1) + eggs
  K <- - (lambda * L) / (lambda - L*molt*surv)

  # write JSON config files
  aqua_pars <- list(
    stochastic = FALSE,
    molt = molt,
    surv = surv,
    K = K,
    L = L
  )

  aqua_pars_path <- tempfile(pattern = "aqua_BH", fileext = ".json")
  jsonlite::write_json(x = aqua_pars, path = aqua_pars_path, pretty = TRUE, digits = 8)

  adult_pars <- list(
    stochastic = FALSE,
    f = f,
    q = q,
    eip = 10,
    p = p,
    psi = diag(1),
    nu = nu,
    M = M,
    Y = 0,
    Z = 0
  )

  adult_pars_path <- tempfile(pattern = "adult_RM", fileext = ".json")
  jsonlite::write_json(x = adult_pars, path = adult_pars_path, pretty = TRUE, digits = 8)

  par <- list(
    p = patches,
    tmax = tmax,
    aqua_path = aqua_pars_path,
    aqua_model = "BH",
    adult_path = adult_pars_path,
    adult_model = "RM"
  )

  global_pars_path <- tempfile(pattern = "config", fileext = ".json")
  jsonlite::write_json(x = par, path = global_pars_path, digits = 8, pretty = TRUE)

  # test API
  root_path <- "http://localhost"

  api <- callr::r_bg(
    function() {
      pr <- plumber::plumb(file = system.file("plumber", "mosquito", "plumber.R", package = "MicroMoB"))
      pr$setDocs(FALSE)
      pr$run(port = 8000)
    }
  )

  Sys.sleep(1)

  expect_true(api$is_alive())

  # test things error properly
  r <- httr::GET(root_path, port = 8000, path = "step_aqua")
  expect_equal(r$status_code, 400)
  expect_equal(httr::content(r), list(error = "model object not found"))

  r <- httr::GET(root_path, port = 8000, path = "step_adult")
  expect_equal(r$status_code, 400)
  expect_equal(httr::content(r), list(error = "model object not found"))

  r <- httr::GET(root_path, port = 8000, path = "config_model_object_mosquito")
  expect_equal(r$status_code, 400)
  expect_equal(httr::content(r), list(error = "model parameters not yet specified"))

  r <- httr::GET(root_path, port = 8000, path = "parameters_adults")
  expect_equal(r$status_code, 400)
  expect_equal(httr::content(r), list(error = "model parameters not yet specified"))

  r <- httr::GET(root_path, port = 8000, path = "parameters_aqua")
  expect_equal(r$status_code, 400)
  expect_equal(httr::content(r), list(error = "model parameters not yet specified"))

  r <- httr::GET(root_path, port = 8000, path = "step_adult")
  expect_equal(r$status_code, 400)
  expect_equal(httr::content(r), list(error = "model object not found"))

  r <- httr::GET(root_path, port = 8000, path = "step_aqua")
  expect_equal(r$status_code, 400)
  expect_equal(httr::content(r), list(error = "model object not found"))

  r <- httr::GET(root_path, port = 8000, path = "output_adult")
  expect_equal(r$status_code, 500)

  r <- httr::GET(root_path, port = 8000, path = "output_aqua")
  expect_equal(r$status_code, 500)

  r <- httr::GET(root_path, port = 8000, path = "clock_tick")
  expect_equal(r$status_code, 400)
  expect_equal(httr::content(r), list(error = "model object not found"))

  # set up model
  r <- httr::GET(root_path, port = 8000, path = "config_mosquito", query = list(path = global_pars_path))
  expect_equal(r$status_code, 200)
  expect_equal(httr::content(r), list(msg = "model parameters successfully read in"))

  r <- httr::GET(root_path, port = 8000, path = "config_model_object_mosquito")
  expect_equal(r$status_code, 200)
  expect_equal(httr::content(r), list(msg = "model successfully set up"))

  r <- httr::GET(root_path, port = 8000, path = "parameters_adults")
  expect_equal(r$status_code, 200)
  expect_true(length(httr::content(r)) == length(adult_pars))

  r <- httr::GET(root_path, port = 8000, path = "parameters_aqua")
  expect_equal(r$status_code, 200)
  expect_true(length(httr::content(r)) == length(aqua_pars))

  for (t in 1:tmax) {
    r <- httr::GET(root_path, port = 8000, path = "step_aqua")
    expect_equal(r$status_code, 200)
    r <- httr::GET(root_path, port = 8000, path = "step_adult")
    expect_equal(r$status_code, 200)
    r <- httr::GET(root_path, port = 8000, path = "clock_tick")
    expect_equal(r$status_code, 200)
  }

  r <- httr::GET(root_path, port = 8000, path = "output_aqua")
  L_out <- httr::content(r)
  r <- httr::GET(root_path, port = 8000, path = "output_adult")
  A_out <- httr::content(r)

  L_out <- as.data.frame(L_out)
  A_out <- as.data.frame(A_out)

  expect_equal(L_out$L, L)
  expect_equal(A_out$M, M)
  expect_equal(A_out$Y, 0)
  expect_equal(A_out$Z, 0)

  api$kill(close_connections = TRUE)
  unlink(x = global_pars_path)
  unlink(x = aqua_pars_path)
  unlink(x = adult_pars_path)

})
