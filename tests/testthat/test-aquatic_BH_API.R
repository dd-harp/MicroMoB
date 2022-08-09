test_that("test aqua BH parameter API", {
  p <- 3
  tmax <- 20
  molt <-  0.5
  surv <- 0.5
  K <-  100
  L <- c(800, 1250, 2500)
  eggs <- c(10, 10, 10)

  mod <- make_MicroMoB(tmax = tmax, p = p)
  setup_aqua_BH(model = mod, stochastic = FALSE, molt = molt, surv = surv, K = K, L = L)

  # change molt
  set_molt_aqua_BH(model = mod, molt = 0.2, times = 2, places = 1)
  expect_equal(get_molt_aqua_BH(model = mod, times = 2, places = 1), 0.2)
  expect_true(all(as.vector(mod$aqua$molt)[-4] == 0.5))

  set_molt_aqua_BH(model = mod, molt = 0.4, times = 5:6, places = 2:3)
  expect_true(all(get_molt_aqua_BH(model = mod, times = 5:6, places = 2:3) == 0.4))
  expect_equal(sum(mod$aqua$molt == 0.5), 55)

  set_molt_aqua_BH(model = mod, molt = matrix(c(0.1,0.2,0.3,0.4),2,2), times = c(10, 20), places = c(1, 3))
  expect_true(all(get_molt_aqua_BH(model = mod, times = c(10,20),places = c(1,3)) == c(0.1,0.2,0.3,0.4)))
  expect_equal(sum(mod$aqua$molt == 0.5), 51)

  expect_error(set_molt_aqua_BH(model = mod, molt = -0.5, times = 1, places = 1))
  expect_error(set_molt_aqua_BH(model = mod, molt = 0.5, times = -5, places = 1))
  expect_error(set_molt_aqua_BH(model = mod, molt = 0.5, times = 5, places = -1))
  expect_error(set_molt_aqua_BH(model = mod, molt = 0.5, times = 5, places = 100))
  expect_error(set_molt_aqua_BH(model = mod, molt = 0.5, times = 100, places = 1))

  # change surv
  set_surv_aqua_BH(model = mod, surv = 0.2, times = 2, places = 1)
  expect_equal(get_surv_aqua_BH(model = mod, times = 2, places = 1), 0.2)
  expect_true(all(as.vector(mod$aqua$surv)[-4] == 0.5))

  set_surv_aqua_BH(model = mod, surv = 0.4, times = 5:6, places = 2:3)
  expect_true(all(get_surv_aqua_BH(model = mod, times = 5:6, places = 2:3) == 0.4))
  expect_equal(sum(mod$aqua$surv == 0.5), 55)

  set_surv_aqua_BH(model = mod, surv = matrix(c(0.1,0.2,0.3,0.4),2,2), times = c(10, 20), places = c(1, 3))
  expect_true(all(get_surv_aqua_BH(model = mod, times = c(10,20),places = c(1,3)) == c(0.1,0.2,0.3,0.4)))
  expect_equal(sum(mod$aqua$surv == 0.5), 51)

  expect_error(set_surv_aqua_BH(model = mod, surv = -0.5, times = 1, places = 1))
  expect_error(set_surv_aqua_BH(model = mod, surv = 0.5, times = -5, places = 1))
  expect_error(set_surv_aqua_BH(model = mod, surv = 0.5, times = 5, places = -1))
  expect_error(set_surv_aqua_BH(model = mod, surv = 0.5, times = 5, places = 100))
  expect_error(set_surv_aqua_BH(model = mod, surv = 0.5, times = 100, places = 1))

  # change K
  set_K_aqua_BH(model = mod, K = 0.2, times = 2, places = 1)
  expect_equal(get_K_aqua_BH(model = mod, times = 2, places = 1), 0.2)
  expect_true(all(as.vector(mod$aqua$K)[-4] == 100))

  set_K_aqua_BH(model = mod, K = 0.4, times = 5:6, places = 2:3)
  expect_true(all(get_K_aqua_BH(model = mod, times = 5:6, places = 2:3) == 0.4))
  expect_equal(sum(mod$aqua$K == 100), 55)

  set_K_aqua_BH(model = mod, K = matrix(c(0.1,0.2,0.3,0.4),2,2), times = c(10, 20), places = c(1, 3))
  expect_true(all(get_K_aqua_BH(model = mod, times = c(10,20),places = c(1,3)) == c(0.1,0.2,0.3,0.4)))
  expect_equal(sum(mod$aqua$K == 100), 51)

  expect_error(set_K_aqua_BH(model = mod, K = -0.5, times = 1, places = 1))
  expect_error(set_K_aqua_BH(model = mod, K = 0.5, times = -5, places = 1))
  expect_error(set_K_aqua_BH(model = mod, K = 0.5, times = 5, places = -1))
  expect_error(set_K_aqua_BH(model = mod, K = 0.5, times = 5, places = 100))
  expect_error(set_K_aqua_BH(model = mod, K = 0.5, times = 100, places = 1))
})
