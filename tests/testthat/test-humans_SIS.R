test_that("human object setup is working", {

  n <- 3
  p <- 3
  tmax <- 10

  mod <- make_MicroMoB(tmax = tmax, p = p)

  theta <- matrix(rexp(9), n, p)
  theta <- theta / rowSums(theta)
  wf <- rep(1, n)
  H <- c(100, 80, 50)
  X <- c(20, 15, 5)
  c <- 0.15

  setup_humans_SIS(model = mod, stochastic = FALSE, theta = theta, H = H, X = X)
  expect_equal(compute_x(mod), (X / H) * c)
  expect_equal(compute_wf(mod), wf)
  expect_equal(compute_Psi(mod), theta)
  expect_equal(compute_H(mod), H)

})

test_that("computation of human biting terms is working", {

  n_strata <- 3
  p <- 2
  n <- p * n_strata
  tmax <- 10

  mod <- make_MicroMoB(tmax = tmax, p = p)

  theta <- matrix(rexp(n*p), n, p)
  theta <- theta / rowSums(theta)
  wf <- rexp(n)
  H <- rpois(n = n, lambda = 100)
  X <- rpois(n = n, lambda = 15)
  c <- 0.15

  setup_humans_SIS(model = mod, stochastic = FALSE, theta = theta, H = H, X = X, wf = wf)

  beta_manual <- matrix(NaN, nrow = n, ncol = p)
  for (i in 1:n) {
    for (j in 1:p) {
      beta_manual[i,j] = (theta[i,j]*wf[i]) / sum(wf * H * theta[, j])
    }
  }

  wf <- compute_wf(mod)
  Psi <- compute_Psi(mod)
  W <- as.vector(t(Psi) %*% (wf * H))
  beta <- diag(wf) %*% Psi %*% diag(1/W)

  expect_equal(beta, beta_manual)

})


test_that("deterministic updates of human SIS model work", {

  tmax <- 1e2
  n <- 3
  p <- 3

  mod <- make_MicroMoB(tmax = tmax, p = p)

  b <- 0.55
  c <- 0.15
  r <- 1/5

  theta <- matrix(rexp(9), n, p)
  theta <- theta / rowSums(theta)
  wf <- rep(1, 3)
  H <- c(100, 80, 50)
  X <- c(20, 15, 5)

  setup_humans_SIS(model = mod, stochastic = FALSE, theta = theta, H = H, X = X, b = b, c = c, r = r)

  for (i in 1:tmax) {
    step_humans(model = mod)
    mod$global$tnow <- mod$global$tnow + 1L
  }

  expect_true(all(mod$human$X < 1e-6))

})


test_that("deterministic updates of human SIS model work with pulsed h", {

  tmax <- 1e2
  n <- 3
  p <- 3

  mod <- make_MicroMoB(tmax = tmax, p = p)

  b <- 1
  c <- 0.15
  r <- 1/200

  theta <- matrix(rexp(9), n, p)
  theta <- theta / rowSums(theta)
  wf <- rep(1, 3)
  H <- c(100, 80, 50)
  X <- c(20, 15, 5)

  setup_humans_SIS(model = mod, stochastic = FALSE, theta = theta, H = H, X = X, b = b, c = c, r = r)

  mod$human$EIR <- rep(qexp(p = 0.25), n)
  step_humans(model = mod)

  expect_true(all(mod$human$X > X))

  mod$global$tnow <- mod$global$tnow + 1L
  prev_X <- mod$human$X

  pulse_X <- mod$human$X

  mod$human$EIR <- rep(0, n)
  for (i in 2:tmax) {
    step_humans(model = mod)
    expect_true(all(mod$human$X < prev_X))
    prev_X <- mod$human$X
  }

  pulse_X_decay <- pulse_X * (1 - pexp(q = r))^(tmax-1)
  expect_equal(mod$human$X, pulse_X_decay)

})


test_that("stochastic updates of human SIS model work", {

  tmax <- 1e2
  n <- 3
  p <- 3

  mod <- make_MicroMoB(tmax = tmax, p = p)

  b <- 0.55
  c <- 0.15
  r <- 1/5

  theta <- matrix(rexp(9), n, p)
  theta <- theta / rowSums(theta)
  wf <- rep(1, 3)
  H <- c(100, 80, 50)
  X <- c(20, 15, 5)

  setup_humans_SIS(model = mod, stochastic = TRUE, theta = theta, H = H, X = X, b = b, c = c, r = r)

  for (i in 1:tmax) {
    step_humans(model = mod)
    mod$global$tnow <- mod$global$tnow + 1L
  }

  expect_true(all(mod$human$X == 0))

})


test_that("stochastic updates of human SIS model work with pulsed h", {

  tmax <- 1e2
  n <- 3
  p <- 3

  mod <- make_MicroMoB(tmax = tmax, p = p)

  b <- 1
  c <- 0.15
  r <- 1/200

  theta <- matrix(rexp(9), n, p)
  theta <- theta / rowSums(theta)
  wf <- rep(1, 3)
  H <- c(100, 80, 50)
  X <- c(20, 15, 5)

  setup_humans_SIS(model = mod, stochastic = TRUE, theta = theta, H = H, X = X, b = b, c = c, r = r)

  mod$human$EIR <- rep(qexp(p = 0.25), n)
  step_humans(model = mod)

  expect_true(all(mod$human$X >= X))

  mod$global$tnow <- mod$global$tnow + 1L
  prev_X <- mod$human$X

  pulse_X <- mod$human$X

  mod$human$EIR <- rep(0, n)
  for (i in 2:tmax) {
    step_humans(model = mod)
    expect_true(all(mod$human$X <= prev_X))
    prev_X <- mod$human$X
  }

})


test_that("test JSON config working", {

  library(jsonlite)

  n <- 6 # number of human population strata
  p <- 2 # number of patches

  theta <- matrix(rexp(n*p), nrow = n, ncol = p)
  theta <- theta / rowSums(theta)
  H <- rep(10, n)
  X <- rep(3, n)

  # sending to JSON does not change R type when read back in
  par <- list(
    "stochastic" = FALSE,
    "theta" = theta,
    "wf" = rep(1, n),
    "H" = H,
    "X" = X,
    "b" = 0.55,
    "c" = 0.15,
    "r" = 1/200
  )

  json_path <- tempfile(pattern = "human_par", fileext = ".json")
  write_json(x = par, path = json_path, digits = NA)
  par_in <- get_config_humans_SIS(path = json_path)
  expect_true(all.equal(par, par_in))

  # reject obviously bad input
  par <- list(
    "stochastic" = "FALSE",
    "theta" = theta,
    "wf" = rep(1, n),
    "H" = H,
    "X" = X,
    "b" = 0.55,
    "c" = 0.15,
    "r" = 1/200
  )

  json_path <- tempfile(pattern = "aqua_par", fileext = ".json")
  write_json(x = par, path = json_path, digits = NA)
  expect_error(get_config_humans_SIS(path = json_path))

  unlink(x = json_path)

})

