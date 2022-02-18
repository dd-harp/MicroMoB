test_that("deterministic MOI model updates correctly, MOI matrix not completely filled", {
  n <- 2

  MOI <- matrix(data = rpois(n = 5*n, lambda = 1e3), nrow = 5, ncol = n)
  MOI <- rbind(MOI, 0)

  EIR <- runif(n = 2, min = 0, max = 0.25)
  b <- 0.55
  h <- EIR * b
  h <- pexp(q = h)

  r <- 1/200
  sigma <- 1
  rho <- r * 1:(nrow(MOI)-1)^sigma
  rho <- pexp(q = rho)

  new_infections <- MOI %*% diag(h)
  recoveries <- diag(rho) %*% MOI[-1, ]

  new_MOI <- MOI
  new_MOI <- new_MOI - new_infections + rbind(0, new_infections[-6, ]) - rbind(0, recoveries) + rbind(recoveries, 0)

  expect_equal(colSums(new_MOI), colSums(MOI))

  mod <- make_MicroMoB(tmax = 1, p = 1)
  setup_humans_MOI(model = mod, stochastic = FALSE, theta = matrix(1, nrow = 2, ncol = 1), H = colSums(MOI), MOI = MOI, b = b, r = r, sigma = sigma)

  mod$human$EIR <- EIR
  step_humans(model = mod)

  expect_equal(mod$human$MOI, new_MOI)

})


test_that("deterministic MOI model updates correctly, MOI matrix is filled", {
  n <- 2

  MOI <- matrix(data = rpois(n = 5*n, lambda = 1e3), nrow = 5, ncol = n)

  EIR <- runif(n = 2, min = 0, max = 0.25)
  b <- 0.55
  h <- EIR * b
  h <- pexp(q = h)

  r <- 1/200
  sigma <- 1
  rho <- r * 1:(nrow(MOI)-1)^sigma
  rho <- pexp(q = rho)

  new_infections <- MOI %*% diag(h)
  recoveries <- diag(rho) %*% MOI[-1, ]

  new_MOI <- MOI
  new_MOI <- rbind(new_MOI, 0)
  new_infections <- rbind(new_infections, 0)
  recoveries <- rbind(recoveries, 0)
  new_MOI <- new_MOI - new_infections + rbind(0, new_infections[-6, ]) - rbind(0, recoveries) + rbind(recoveries, 0)

  expect_equal(colSums(new_MOI), colSums(MOI))

  mod <- make_MicroMoB(tmax = 1, p = 1)
  setup_humans_MOI(model = mod, stochastic = FALSE, theta = matrix(1, nrow = 2, ncol = 1), H = colSums(MOI), MOI = MOI, b = b, r = r, sigma = sigma)

  mod$human$EIR <- EIR
  step_humans(model = mod)

  expect_equal(mod$human$MOI, new_MOI)

})


test_that("stochastic MOI model updates correctly, MOI matrix not completely filled", {
  n <- 2

  MOI <- matrix(data = rpois(n = 5*n, lambda = 1e6), nrow = 5, ncol = n)
  MOI <- rbind(MOI, 0)

  EIR <- runif(n = 2, min = 0, max = 0.25)
  r <- 1/50

  mod <- make_MicroMoB(tmax = 1, p = 1)
  setup_humans_MOI(model = mod, stochastic = FALSE, theta = matrix(1, nrow = 2, ncol = 1), H = colSums(MOI), MOI = MOI, r = r)

  mod$human$EIR <- EIR
  step_humans(model = mod)
  MOI_det <- mod$human$MOI

  mod <- make_MicroMoB(tmax = 1, p = 1)
  setup_humans_MOI(model = mod, stochastic = TRUE, theta = matrix(1, nrow = 2, ncol = 1), H = colSums(MOI), MOI = MOI, r = r)

  mod$human$EIR <- EIR
  step_humans(model = mod)
  MOI_sto <- mod$human$MOI

  diffs <- abs(MOI_det - MOI_sto) / MOI_det
  expect_true(all(diffs < 0.10))

  # pval <- ks.test(as.vector(MOI_det), as.vector(MOI_sto))$p.value
  # expect_true(pval > 0.8)
})


test_that("stochastic MOI model updates correctly, MOI matrix is filled", {
  n <- 2

  MOI <- matrix(data = rpois(n = 5*n, lambda = 1e6), nrow = 5, ncol = n)

  EIR <- runif(n = 2, min = 0, max = 0.25)
  r <- 1/50

  mod <- make_MicroMoB(tmax = 1, p = 1)
  setup_humans_MOI(model = mod, stochastic = FALSE, theta = matrix(1, nrow = 2, ncol = 1), H = colSums(MOI), MOI = MOI, r = r)

  mod$human$EIR <- EIR
  step_humans(model = mod)
  MOI_det <- mod$human$MOI

  mod <- make_MicroMoB(tmax = 1, p = 1)
  setup_humans_MOI(model = mod, stochastic = TRUE, theta = matrix(1, nrow = 2, ncol = 1), H = colSums(MOI), MOI = MOI, r = r)

  mod$human$EIR <- EIR
  step_humans(model = mod)
  MOI_sto <- mod$human$MOI

  diffs <- abs(MOI_det - MOI_sto) / MOI_det
  expect_true(all(diffs < 0.10))

  # pval <- ks.test(as.vector(MOI_det), as.vector(MOI_sto))$p.value
  # expect_true(pval > 0.8)
})


test_that("test JSON config working", {

  library(jsonlite)

  n <- 6 # number of human population strata
  p <- 2 # number of patches

  theta <- matrix(rexp(n*p), nrow = n, ncol = p)
  theta <- theta / rowSums(theta)
  H <- rep(10, n)
  MOI <- matrix(0, nrow = 10, ncol = n)
  MOI[1, ] <- H

  # sending to JSON does not change R type when read back in
  par <- list(
    "stochastic" = FALSE,
    "theta" = theta,
    "wf" = rep(1, n),
    "H" = H,
    "MOI" = MOI,
    "b" = 0.55,
    "c" = 0.15,
    "r" = 1/200,
    "sigma" = 1
  )

  json_path <- tempfile(pattern = "human_par", fileext = ".json")
  write_json(x = par, path = json_path, digits = NA)
  par_in <- get_config_humans_MOI(path = json_path)
  expect_true(all.equal(par, par_in))

  # reject obviously bad input
  par <- list(
    "stochastic" = FALSE,
    "theta" = theta[1],
    "wf" = rep(1, n),
    "H" = H,
    "MOI" = MOI,
    "b" = 0.55,
    "c" = 0.15,
    "r" = 1/200,
    "sigma" = 1
  )

  json_path <- tempfile(pattern = "aqua_par", fileext = ".json")
  write_json(x = par, path = json_path, digits = NA)
  expect_error(get_config_humans_MOI(path = json_path))

  unlink(x = json_path)

})


test_that("JSON parameters can read in", {
  path <- system.file("extdata", "humans_MOI.json", package = "MicroMoB")
  pars <- get_config_humans_MOI(path = path)

  expect_true(length(pars) == 9L)

  expect_true(is.logical(pars$stochastic))
  stopifnot(length(pars$stochastic) == 1L)

  expect_true(is.numeric(pars$theta))
  expect_true(is.matrix(pars$theta))

  expect_true(is.numeric(pars$wf))
  expect_true(is.vector(pars$wf))

  expect_true(is.numeric(pars$H))
  expect_true(is.vector(pars$H))

  expect_true(is.numeric(pars$MOI))
  expect_true(is.matrix(pars$MOI))

  expect_true(is.numeric(pars$b))
  expect_true(is.vector(pars$b))

  expect_true(is.numeric(pars$c))
  expect_true(is.vector(pars$c))

  expect_true(is.numeric(pars$r))
  expect_true(is.vector(pars$r))

  expect_true(is.numeric(pars$sigma))
  expect_true(is.vector(pars$sigma))
})
