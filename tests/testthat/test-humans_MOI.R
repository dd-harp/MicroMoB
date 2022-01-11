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


