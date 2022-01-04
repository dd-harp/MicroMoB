test_that("MOI model updates correctly, MOI matrix not completely filled", {
  n <- 2

  MOI <- matrix(data = rpois(n = 5*n, lambda = 1e3), nrow = 5, ncol = n)
  MOI <- rbind(MOI, 0)

  EIR <- runif(n = 2, min = 0, max = 0.25)
  b <- 0.55
  h <- EIR * b

  r <- 1/200
  sigma <- 1
  rho <- r * 1:(nrow(MOI)-1)^sigma

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
