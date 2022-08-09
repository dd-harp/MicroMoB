test_that("test mosquito RM parameter API", {
  tmax <- 20

  f <- 0.3
  q <- 1
  a <- f * q
  psi <- matrix(
    c(0.9, 0.025, 0.075,
      0.15, 0.6, 0.25,
      0.01, 0.04, 0.95), nrow = 3, ncol = 3, byrow = TRUE
  )
  M <- c(100, 150, 120)
  Y <- c(10, 15, 8)
  Z <- c(8, 10, 4)

  mod <- make_MicroMoB(tmax = tmax, p = 3)
  setup_mosquito_RM(mod, stochastic = FALSE, f = f, q = q, eip = 5, p = 0.9, psi = psi, M = M, Y = Y, Z = Z)

  # change f
  expect_error(modify_f_mosquito_RM(model = mod, f = -1))
  expect_error(modify_f_mosquito_RM(model = mod, f = rep(1, 100)))
  expect_error(modify_f_mosquito_RM(model = mod, f = '5'))
  expect_error(modify_f_mosquito_RM(model = mod, f = Inf))
  expect_error(modify_f_mosquito_RM(model = mod, f = NA))
  expect_error(modify_f_mosquito_RM(model = mod, f = NaN))
  expect_error(modify_f_mosquito_RM(model = mod, f = NULL))
  modify_f_mosquito_RM(model = mod, f = 1)
  expect_equal(mod$mosquito$f, 1)

  # change q
  expect_error(modify_q_mosquito_RM(model = mod, q = -1))
  expect_error(modify_q_mosquito_RM(model = mod, q = rep(1, 100)))
  expect_error(modify_q_mosquito_RM(model = mod, q = '5'))
  expect_error(modify_q_mosquito_RM(model = mod, q = Inf))
  expect_error(modify_q_mosquito_RM(model = mod, q = NA))
  expect_error(modify_q_mosquito_RM(model = mod, q = NaN))
  expect_error(modify_q_mosquito_RM(model = mod, q = NULL))
  modify_q_mosquito_RM(model = mod, q = 1)
  expect_equal(mod$mosquito$q, 1)

  # change EIP
  expect_error(modify_eip_mosquito_RM(model = mod, eip = 10, times = 100))
  expect_error(modify_eip_mosquito_RM(model = mod, eip = 10, times = 0))
  expect_error(modify_eip_mosquito_RM(model = mod, eip = -10, times = 5))
  expect_error(modify_eip_mosquito_RM(model = mod, eip = NaN, times = 5))
  expect_error(modify_eip_mosquito_RM(model = mod, eip = c(10, 50, 2), times = c(1, 5)))

  modify_eip_mosquito_RM(model = mod, eip = 10, times = c(1,5))
  expect_true(all(mod$mosquito$eip[c(1,5)] == 10))
  expect_true(all(mod$mosquito$eip[!c(1,5)] == 5))
  expect_true(all(dim(mod$mosquito$ZZ_shift) == c(10, 10)))
  expect_true(nrow(mod$mosquito$ZZ) == 10)
  expect_true(mod$mosquito$maxEIP == 10)

  modify_eip_mosquito_RM(model = mod, eip = c(10, 20), times = c(1,5))
  expect_true(all(mod$mosquito$eip[c(1,5)] == c(10, 20)))
  expect_true(all(mod$mosquito$eip[!c(1,5)] == 5))
  expect_true(all(dim(mod$mosquito$ZZ_shift) == c(20, 20)))
  expect_true(nrow(mod$mosquito$ZZ) == 20)
  expect_true(mod$mosquito$maxEIP == 20)

  # change p
  modify_p_mosquito_RM(model = mod, p = 0.5, times = 2, places = 1)
  expect_equal(mod$mosquito$p[1, 2], 0.5)
  expect_true(all(as.vector(mod$mosquito$p)[-4] == 0.9))

  modify_p_mosquito_RM(model = mod, p = 0.5, times = 5:6, places = 2:3)
  expect_true(all(mod$mosquito$p[2:3, 5:6] == 0.5))
  expect_equal(sum(mod$mosquito$p == 0.9), 55)

  modify_p_mosquito_RM(model = mod, p = matrix(c(0.1,0.2,0.3,0.4),2,2), times = c(10, 20), places = c(1, 3))
  expect_true(all(mod$mosquito$p[c(1,3),c(10,20)] == c(0.1,0.2,0.3,0.4)))
  expect_equal(sum(mod$mosquito$p == 0.9), 51)

  expect_error(modify_p_mosquito_RM(model = mod, p = -0.5, times = 1, places = 1))
  expect_error(modify_p_mosquito_RM(model = mod, p = 0.5, times = -5, places = 1))
  expect_error(modify_p_mosquito_RM(model = mod, p = 0.5, times = 5, places = -1))
  expect_error(modify_p_mosquito_RM(model = mod, p = 0.5, times = 5, places = 100))
  expect_error(modify_p_mosquito_RM(model = mod, p = 0.5, times = 100, places = 1))

  # change psi
  expect_error(modify_psi_mosquito_RM(model = mod, psi = 1))
  expect_error(modify_psi_mosquito_RM(model = mod, psi = psi * -1))
  expect_error(modify_psi_mosquito_RM(model = mod, psi = psi * 5))

  modify_psi_mosquito_RM(model = mod, psi = diag(3))
  expect_equal(mod$mosquito$psi, diag(3))

  # change nu
  expect_error(modify_nu_mosquito_RM(model = mod, nu = -1))
  expect_error(modify_nu_mosquito_RM(model = mod, nu = rep(1, 100)))
  expect_error(modify_nu_mosquito_RM(model = mod, nu = '5'))
  expect_error(modify_nu_mosquito_RM(model = mod, nu = Inf))
  expect_error(modify_nu_mosquito_RM(model = mod, nu = NA))
  expect_error(modify_nu_mosquito_RM(model = mod, nu = NaN))
  expect_error(modify_nu_mosquito_RM(model = mod, nu = NULL))
  modify_nu_mosquito_RM(model = mod, nu = 1)
  expect_equal(mod$mosquito$nu, 1)

  # change kappa
  expect_error(modify_kappa_mosquito_RM(model = mod, kappa = -1))
  expect_error(modify_kappa_mosquito_RM(model = mod, kappa = rep(1, 100)))
  expect_error(modify_kappa_mosquito_RM(model = mod, kappa = '5'))
  expect_error(modify_kappa_mosquito_RM(model = mod, kappa = Inf))
  expect_error(modify_kappa_mosquito_RM(model = mod, kappa = NA))
  expect_error(modify_kappa_mosquito_RM(model = mod, kappa = NaN))
  expect_error(modify_kappa_mosquito_RM(model = mod, kappa = NULL))
  modify_kappa_mosquito_RM(model = mod, kappa = c(0.5, 0.6, 0.3))
  expect_equal(mod$mosquito$kappa, c(0.5, 0.6, 0.3))
})
