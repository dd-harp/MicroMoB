test_that("test bloodmeal with simple RM setup", {
  patches <- 1
  n <- 1
  tmax <- 5

  # human parameters
  theta <- diag(n)
  H <- 100
  X <- 30
  b <- 0.55
  c <- 0.15
  r <- 1/200

  # mosquito parameters
  f <- 0.3
  q <- 1
  eip <- 10
  p <- 0.9
  M <- 500
  Y <- 100
  Z <- 80
  psi <- diag(patches)

  mod <- make_MicroMoB(tmax = tmax, p = patches)
  setup_humans_SIS(mod, stochastic = FALSE, theta = theta, H = H, X = X, b = b, c = c, r = r)
  setup_aqua_trace(mod, stochastic = FALSE, lambda = 0)
  setup_mosquito_RM(mod, stochastic = FALSE, f = f, q = q, eip = eip, p = p, psi = psi, M = M, Y = Y, Z = Z)
  setup_alternative_trace(mod)
  setup_visitor_trace(mod)

  # compute EIR by hand
  EIR <- 1/H * f * q * Z

  # compute kappa by hand
  kappa <- c * (X/H)

  # compute EIR and kappa in MicroMoB
  compute_bloodmeal(mod)

  # check
  expect_equal(as.vector(mod$human$EIR), EIR)
  expect_equal(as.vector(mod$mosquito$kappa), kappa)

})
