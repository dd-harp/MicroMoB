test_that("setting up time spent (day) works", {

  model <- setup_model_object()
  setup_human("strata", model = model, H = rep(10, 3))
  expect_error(setup_timespent("blargh", model = model))

  # no strata, just patches
  model <- setup_model_object()
  setup_human("strata", model = model, H = rep(10, 3))
  setup_timespent("day", model = model)
  setup_biteweight("simple", model = model)

  expect_true(inherits(model$human$timespent, "day"))
  expect_equal(model$human$timespent$Theta_t, diag(3))

  Psi_t <- compute_Psi(human = model$human, xi = 1, t = 1)
  W <- compute_W(human = model$human, Psi_t = Psi_t, t = 1)

  expect_equal(W, model$human$H)

  # 2 patch, 3 strata
  H <- c(50, 20, 10, 30, 20, 10)
  J <- matrix(
    c(1, 0, 1, 0, 1, 0,
      0, 1, 0, 1, 0, 1),
    nrow = 2, ncol = 6, byrow = TRUE
  )

  # theta keeps everyone in one place
  model <- setup_model_object()
  setup_human("strata", model = model, H = H, J = J)
  setup_timespent("day", model = model)

  expect_equal(model$human$timespent$Theta_t %*% H, model$human$J %*% model$human$H)
  expect_equal(sum(model$human$timespent$Theta_t %*% H), sum(H))
})


test_that("setting up time spent (dt) works", {

  # 2 patch, 3 strata
  H <- c(50, 20, 10, 30, 20, 10)
  J <- matrix(
    c(1, 0, 1, 0, 1, 0,
      0, 1, 0, 1, 0, 1),
    nrow = 2, ncol = 6, byrow = TRUE
  )

  expect_equal(as.vector(J %*% H), c(sum(H[1], H[3], H[5]), sum(H[2], H[4], H[6])))

  theta_day <- matrix(
    c(0.4, 0.6,
      0.5, 0.5,
      0.2, 0.8,
      0.5, 0.5,
      0.6, 0.4,
      0.2, 0.8), nrow = length(H), ncol = 2, byrow = TRUE
  )

  theta_night <- matrix(
    c(0.9, 0.1,
      0.6, 0.4,
      0.5, 0.5,
      0.6, 0.4,
      0.8, 0.2,
      0.4, 0.6), nrow = length(H), ncol = 2, byrow = TRUE
  )

  theta <- list(theta_day, theta_night)

  model <- setup_model_object()
  setup_human("strata", model = model, H = H, J = J)
  setup_timespent("dt", model = model, theta = theta)
  setup_biteweight("simple", model = model)

  xi <- c(0.15, 0.85)
  Psi_t <- compute_Psi(human = model$human, xi = xi, t = 1)
  W <- compute_W(human = model$human, Psi_t = Psi_t, t = 1)

  wf <- rep(1, length(H))
  W_expected <- ((t(theta_day) %*% (wf * H)) * (xi[1])) + ((t(theta_night) %*% (wf * H)) * (xi[2]))

  expect_equal(rowSums(W), as.vector(W_expected))

  # 3 strata, 2 patch
  H_strata <- c(70, 40, 30)
  J_strata <- matrix(
    c(50/70,  10/40, 20/30,
      1-50/70, 1-10/40, 1-20/30),
    nrow = 2,
    ncol = 3, byrow = TRUE
  )

  residency <- strata_to_residency_proportion(H_strata = H_strata, J_strata = J_strata)
  J <- residency$J
  H <- residency$H

  model <- setup_model_object()
  setup_human("strata", model = model, H = H, J = J)
  setup_timespent("dt", model = model, theta = theta)
  setup_biteweight("simple", model = model)

  xi <- c(0.15, 0.85)
  Psi_t <- compute_Psi(human = model$human, xi = xi, t = 1)
  W <- compute_W(human = model$human, Psi_t = Psi_t, t = 1)

  expect_equal(rowSums(W), as.vector(W_expected))
})
