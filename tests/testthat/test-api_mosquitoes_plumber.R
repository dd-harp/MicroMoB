# Setup by starting APIs
root_path <- "http://localhost"

# api_port <- 8000
api_port <- plumber:::findPort()

api <- callr::r_bg(
  function(port) {
    pr <- plumber::plumb(system.file("plumber", "mosquito", "plumber.R", package = "MicroMoB"))
    pr$run(port = port)
  }, args = list(port = api_port)
)

Sys.sleep(3)

withr::defer(api$kill())

test_that("API is alive", {
  expect_true(api$is_alive())
})

test_that("endpoints error properly", {

  r <- httr::GET(root_path, port = api_port, path = "config_global_parameters")
  expect_equal(httr::status_code(r), 500)
  expect_equal(httr::content(r), list(error = "500 - Internal server error"))

  r <- httr::GET(root_path, port = api_port, path = "config_aqua_parameters")
  expect_equal(httr::status_code(r), 400)
  expect_equal(httr::content(r), list(error = "model parameters not yet specified"))

  r <- httr::GET(root_path, port = api_port, path = "config_adult_parameters")
  expect_equal(httr::status_code(r), 400)
  expect_equal(httr::content(r), list(error = "model parameters not yet specified"))

  r <- httr::GET(root_path, port = api_port, path = "parameters_global")
  expect_equal(httr::status_code(r), 400)
  expect_equal(httr::content(r), list(error = "model parameters not yet specified"))

  r <- httr::GET(root_path, port = api_port, path = "parameters_adult")
  expect_equal(httr::status_code(r), 400)
  expect_equal(httr::content(r), list(error = "model parameters not yet specified"))

  r <- httr::GET(root_path, port = api_port, path = "parameters_aqua")
  expect_equal(httr::status_code(r), 400)
  expect_equal(httr::content(r), list(error = "model parameters not yet specified"))

  r <- httr::GET(root_path, port = api_port, path = "setup_model")
  expect_equal(httr::status_code(r), 400)
  expect_equal(httr::content(r), list(error = "model parameters not yet specified"))

  r <- httr::GET(root_path, port = api_port, path = "setup_aqua")
  expect_equal(httr::status_code(r), 400)
  expect_equal(httr::content(r), list(error = "model parameters not yet specified"))

  r <- httr::GET(root_path, port = api_port, path = "setup_adult")
  expect_equal(httr::status_code(r), 400)
  expect_equal(httr::content(r), list(error = "model parameters not yet specified"))

  r <- httr::GET(root_path, port = api_port, path = "step_aqua")
  expect_equal(httr::status_code(r), 400)
  expect_equal(httr::content(r), list(error = "model object not found"))

  r <- httr::GET(root_path, port = api_port, path = "step_adult")
  expect_equal(httr::status_code(r), 400)
  expect_equal(httr::content(r), list(error = "model object not found"))

  r <- httr::GET(root_path, port = api_port, path = "output_aqua")
  expect_equal(httr::status_code(r), 500)

  r <- httr::GET(root_path, port = api_port, path = "output_adult")
  expect_equal(httr::status_code(r), 500)

  r <- httr::GET(root_path, port = api_port, path = "clock_tick")
  expect_equal(r$status_code, 400)
  expect_equal(httr::content(r), list(error = "model object not found"))

})


test_that("run simulation via API", {

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

  # set up model
  r <- httr::GET(root_path, port = api_port, path = "config_global_parameters", query = list(path = global_pars_path))
  expect_equal(r$status_code, 200)
  expect_equal(httr::content(r), list(msg = "global parameters successfully read in"))

  r <- httr::GET(root_path, port = api_port, path = "config_aqua_parameters")
  expect_equal(r$status_code, 200)
  expect_equal(httr::content(r), list(msg = "aquatic component parameters successfully read in"))

  r <- httr::GET(root_path, port = api_port, path = "config_adult_parameters")
  expect_equal(r$status_code, 200)
  expect_equal(httr::content(r), list(msg = "adult component parameters successfully read in"))

  r <- httr::GET(root_path, port = api_port, path = "setup_model")
  expect_equal(r$status_code, 200)
  expect_equal(httr::content(r), list(msg = "model object successfully set up"))

  r <- httr::GET(root_path, port = api_port, path = "setup_aqua")
  expect_equal(r$status_code, 200)
  expect_equal(httr::content(r), list(msg = "aquatic component successfully set up"))

  r <- httr::GET(root_path, port = api_port, path = "setup_adult")
  expect_equal(r$status_code, 200)
  expect_equal(httr::content(r), list(msg = "adult component successfully set up"))

  r <- httr::GET(root_path, port = api_port, path = "parameters_adult")
  expect_equal(r$status_code, 200)
  expect_true(length(httr::content(r)) == length(adult_pars))

  r <- httr::GET(root_path, port = api_port, path = "parameters_aqua")
  expect_equal(r$status_code, 200)
  expect_true(length(httr::content(r)) == length(aqua_pars))

  for (t in 1:tmax) {
    r <- httr::GET(root_path, port = api_port, path = "step_aqua")
    expect_equal(r$status_code, 200)
    r <- httr::GET(root_path, port = api_port, path = "step_adult")
    expect_equal(r$status_code, 200)
    r <- httr::GET(root_path, port = api_port, path = "clock_tick")
    expect_equal(r$status_code, 200)
  }

  r <- httr::GET(root_path, port = api_port, path = "output_aqua")
  L_out <- httr::content(r)
  r <- httr::GET(root_path, port = api_port, path = "output_adult")
  A_out <- httr::content(r)

  L_out <- as.data.frame(L_out)
  A_out <- as.data.frame(A_out)

  expect_equal(L_out$L, L)
  expect_equal(A_out$M, M)
  expect_equal(A_out$Y, 0)
  expect_equal(A_out$Z, 0)

  unlink(x = global_pars_path)
  unlink(x = aqua_pars_path)
  unlink(x = adult_pars_path)
})
