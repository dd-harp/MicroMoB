test_that("Model updates correctly with RM mosquitoes, SIS humans, 2 patches", {

  H_strata <- c(1.2e4, 8e4, 2.5e4)
  J_strata <- matrix(data = c(
    0.5, 0.25, 0.25,
    0.5, 0.75, 0.75
  ), nrow = 2, ncol = 3, byrow = TRUE)

  residency <- strata_to_residency_proportion(H_strata = H_strata, J_strata = J_strata)

  H <- residency$H

  patches <- 2
  n <- length(H)
  tmax <- 5

  # human parameters
  theta <- matrix(c(
    0.9, 0.1,
    0.05, 0.95
  ), nrow = 2, ncol = 2, byrow = TRUE)
  theta <- do.call(rbind, replicate(n = 3, expr = theta, simplify = FALSE))

  X <- rbinom(n = n, size = H, prob = 0.25)
  b <- 0.55
  c <- 0.15
  r <- 1/200

  # mosquito parameters
  f <- 0.3
  q <- 1
  eip <- 10
  p <- 0.9
  M <- rep(5e5, patches)
  Y <- rep(1e5, patches)
  Z <- rep(8e3, patches)
  psi <- matrix(
    c(0.95, 0.05,
      0.05, 0.95),
    nrow = 2, ncol = 2, byrow = TRUE
  )

  # setup model
  mod <- make_MicroMoB(tmax = tmax, p = patches)
  setup_humans_SIS(mod, stochastic = FALSE, theta = theta, H = H, X = X, b = b, c = c, r = r)
  setup_aqua_trace(mod, stochastic = FALSE, lambda = rep(0, patches))
  setup_mosquito_RM(mod, stochastic = FALSE, f = f, q = q, eip = eip, p = p, psi = psi, M = M, Y = Y, Z = Z)
  setup_alternative_trace(mod)
  setup_visitor_trace(mod)

  # calculate 1 step changes
  W <- as.vector(t(theta) %*% H)
  beta <- theta %*% diag(1/W)
  EIR <- as.vector(beta %*% (f*q*Z))
  h <- EIR * b

  x <- (X/H) * c
  kappa <- as.vector(t(beta) %*% (x*H))
  h_mosquito <- kappa * f * q

  # human dynamics
  recoveries <- pexp(q = r) * X
  infections_human <- pexp(q = h) * (H - X)

  # mosquito dynamics
  infections_mosy <- (f*q) * kappa * (M - Y)

  # compute update
  compute_bloodmeal(mod)
  step_humans(mod)
  step_mosquitoes(mod)

  # tests
  expect_equal(mod$human$EIR, EIR)
  expect_equal(mod$human$EIR * mod$human$b, h)
  expect_equal(mod$mosquito$kappa, kappa)
  expect_equal(mod$human$X, X - recoveries + infections_human)
  expect_equal(mod$mosquito$Y, as.vector(((Y + infections_mosy) * p) %*% psi))
  expect_equal(mod$mosquito$M, as.vector((M * p) %*% psi))

})
