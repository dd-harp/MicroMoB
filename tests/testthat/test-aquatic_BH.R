test_that("3 patch BH aqua model works", {
  p <- 3
  tmax <- 1
  molt <-  c(0.2, 0.05, 0.1)
  surv <- c(0.9, 0.95, 0.975)
  K <-  c(1e3, 2.5e3, 5e3)
  L <- c(800, 1250, 2500)
  eggs <- c(10, 10, 10)

  # check deterministic state update works
  mod <- make_MicroMoB(tmax = tmax, p = p)
  setup_aqua_BH(model = mod, stochastic = FALSE, molt = molt, surv = surv, K = K, L = L)
  setup_mosquito_trace(model = mod, oviposit = eggs)

  expect_equal(compute_emergents(mod), rep(0, 3))
  Lt <- eggs + (1-molt)*surv*L*(K/(L+K))
  At <- molt*surv*L*(K/(L+K))

  step_aqua(model = mod)

  expect_equal(mod$aqua$L, Lt)
  expect_equal(compute_emergents(mod), At)

  # check stochastic state update works
  mod <- make_MicroMoB(tmax = tmax, p = p)
  setup_aqua_BH(model = mod, stochastic = TRUE, molt = molt, surv = surv, K = K, L = L)
  setup_mosquito_trace(model = mod, oviposit = eggs)

  expect_equal(compute_emergents(mod), rep(0, 3))

  step_aqua(model = mod)

  expect_true(all(mod$aqua$L > 0))
  expect_equal(order(mod$aqua$L), order(Lt))
  expect_equal(order(compute_emergents(mod)), order(At))
})


test_that("BH equilibrium works", {
  p <- 1
  tmax <- 10

  # lambda out and eggs in are known quantities
  lambda <- 10
  eggs <- 100

  # static pars
  molt <-  0.1
  surv <- 0.9

  # solve L
  L <- lambda * ((1/molt) - 1) + eggs
  K <- - (lambda * L) / (lambda - L*molt*surv)

  # model
  mod <- make_MicroMoB(tmax = tmax, p = p)
  setup_aqua_BH(model = mod, stochastic = FALSE, molt = molt, surv = surv, K = K, L = L)
  setup_mosquito_trace(model = mod, oviposit = eggs)

  while (mod$global$tnow <= tmax) {
    step_aqua(model = mod)
    mod$global$tnow <- mod$global$tnow + 1L
  }

  expect_equal(mod$aqua$L, L)
  expect_equal(compute_emergents(mod), lambda)

})
