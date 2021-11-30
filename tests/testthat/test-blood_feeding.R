test_that("daily beta, biting distribution matrix test", {

  p <- 3
  s <- 2
  n <- s*p

  H_count <- matrix(
    data = c(
      10, 15,
      25, 30,
      20, 18
    ), nrow = p, ncol = s, byrow = TRUE
  )

  H <- as.vector(H_count)

  theta <- matrix(
    data = c(
      0.6, 0.3, 0.1,
      0.2, 0.6, 0.2,
      0.5, 0.3, 0.2,
      0.1, 0.6, 0.3,
      0.6, 0.2, 0.2,
      0.1, 0.7, 0.2
    ), nrow = n, ncol = p, byrow = TRUE
  )

  wf <- c(rep(0.7, p), rep(1.4, p))

  # calculate beta manually
  W <- t(theta) %*% (wf * H)
  beta <- diag(wf) %*% theta %*% diag(1/as.vector(W))

  # calculate beta in Micro-MoB
  residency <- strata_to_residency_counts(H_counts = H_count)

  model <- setup_model_object()
  setup_human("strata", model = model, H = residency$H, J = residency$J)
  setup_timespent("day", model = model, theta = theta)
  setup_biteweight("simple", model = model, wf = wf)

  Psi_t <- compute_Psi(human = model$human, xi = 1, t = 1)
  beta_compare <- compute_beta(human = model$human, Psi_t = Psi_t, t = 1)

  # columns of beta_strata will sum to one
  beta_strata <- diag(H) %*% beta

  expect_equal(beta, beta_compare)
  expect_true(all.equal(colSums(beta_strata), rep(1, 3)))

  # test bites are conserved
  bites <- pmax(1, rpois(n = p, lambda = 50))
  HBR <- beta %*% bites
  bites_strata <- beta_strata %*% bites

  expect_equal(bites_strata / H, HBR)

})


test_that("daily v. fractional (2 chunks) beta, expect exact match with same theta", {

  p <- 3
  s <- 2
  n <- s*p

  H_count <- matrix(
    data = c(
      10, 15,
      25, 30,
      20, 18
    ), nrow = p, ncol = s, byrow = TRUE
  )

  H <- as.vector(H_count)

  theta_day <- matrix(
    data = c(
      0.6, 0.3, 0.1,
      0.2, 0.6, 0.2,
      0.5, 0.3, 0.2,
      0.1, 0.6, 0.3,
      0.6, 0.2, 0.2,
      0.1, 0.7, 0.2
    ), nrow = n, ncol = p, byrow = TRUE
  )

  theta_night <- theta_day
  theta <- list(theta_day, theta_night)

  wf <- c(rep(0.9, p), rep(1.1, p))
  xi <- c(0.25, 0.75)

  # calculate fractional beta manually
  W_day <- (t(theta_day) * xi[1]) %*% (wf * H)
  W_night <- (t(theta_night) * xi[2]) %*% (wf * H)

  beta_dt_day <- diag(wf) %*% (theta_day * xi[1]) %*% diag(1/as.vector(W_day))
  beta_dt_night <- diag(wf) %*% (theta_night * xi[2]) %*% diag(1/as.vector(W_night))
  beta_dt_manual <- (beta_dt_day * xi[1]) + (beta_dt_night * xi[2])

  # calculate fractional beta in Micro-MoB
  residency <- strata_to_residency_counts(H_counts = H_count)

  model <- setup_model_object()
  setup_human("strata", model = model, H = residency$H, J = residency$J)
  setup_timespent("dt", model = model, theta = theta)
  setup_biteweight("simple", model = model, wf = wf)

  Psi_t <- compute_Psi(human = model$human, xi = xi, t = 1)
  beta_dt <- compute_beta(human = model$human, Psi_t = Psi_t, t = 1)

  expect_equal(beta_dt[, , 1], beta_dt_day)
  expect_equal(beta_dt[, , 2], beta_dt_night)

  # calculate daily beta
  model <- setup_model_object()
  setup_human("strata", model = model, H = residency$H, J = residency$J)
  setup_timespent("day", model = model, theta = theta[[1]])
  setup_biteweight("simple", model = model, wf = wf)

  Psi_t <- compute_Psi(human = model$human, xi = 1, t = 1)
  beta_daily <- compute_beta(human = model$human, Psi_t = Psi_t, t = 1)

  # calculate fractional  beta manually
  theta_daily <- theta_day + theta_night
  theta_daily <- diag(1/rowSums(theta_daily)) %*% theta_daily

  W <- t(theta_daily) %*% (wf * H)
  beta_daily_manual <- diag(wf) %*% theta_daily %*% diag(1/as.vector(W))

  # expect daily and fractional to be equal with same theta
  expect_equal(beta_dt_manual, beta_daily)
  expect_equal(beta_daily_manual, beta_daily)

  # test all beta matrices fulfill requirement on columns
  expect_true(all.equal(colSums(diag(H) %*% beta_dt_day), rep(1, 3)))
  expect_true(all.equal(colSums(diag(H) %*% beta_dt_night), rep(1, 3)))
  expect_true(all.equal(colSums(diag(H) %*% beta_dt_manual), rep(1, 3)))
  expect_true(all.equal(colSums(diag(H) %*% beta_dt[, , 1]), rep(1, 3)))
  expect_true(all.equal(colSums(diag(H) %*% beta_dt[, , 2]), rep(1, 3)))
  expect_true(all.equal(colSums(diag(H) %*% beta_daily), rep(1, 3)))
  expect_true(all.equal(colSums(diag(H) %*% beta_daily_manual), rep(1, 3)))

  # bites are conserved
  bites <- pmax(1, rpois(n = p, lambda = 50))

  HBR_day <- beta_dt_day %*% bites
  bites_strata_day <- diag(H) %*% beta_dt_day %*% bites

  expect_equal(bites_strata_day / H, HBR_day)
  expect_equal((diag(H) %*% beta_dt[, , 1] %*% bites) / H, HBR_day)

  HBR_night <- beta_dt_night %*% bites
  bites_strata_night <- diag(H) %*% beta_dt_night %*% bites

  expect_equal(bites_strata_night / H, HBR_night)
  expect_equal((diag(H) %*% beta_dt[, , 2] %*% bites) / H, HBR_night)

})




test_that("daily v. fractional (2 chunks) beta, thetas differ", {

  p <- 3
  s <- 2
  n <- s*p

  H_count <- matrix(
    data = c(
      10, 15,
      25, 30,
      20, 18
    ), nrow = p, ncol = s, byrow = TRUE
  )

  H <- as.vector(H_count)

  theta_day <- matrix(
    data = c(
      0.6, 0.3, 0.1,
      0.2, 0.6, 0.2,
      0.5, 0.3, 0.2,
      0.1, 0.6, 0.3,
      0.6, 0.2, 0.2,
      0.1, 0.7, 0.2
    ), nrow = n, ncol = p, byrow = TRUE
  )

  theta_night <- matrix(
    data = c(
      0.5, 0.3, 0.2,
      0.1, 0.6, 0.3,
      0.4, 0.3, 0.3,
      0.2, 0.5, 0.3,
      0.5, 0.3, 0.2,
      0.1, 0.6, 0.3
    ), nrow = n, ncol = p, byrow = TRUE
  )

  theta <- list(theta_day, theta_night)

  wf <- c(rep(0.9, p), rep(1.1, p))
  xi <- c(1/3, 2/3)

  # calculate fractional beta manually
  W_day <- (t(theta_day) * xi[1]) %*% (wf * H)
  W_night <- (t(theta_night) * xi[2]) %*% (wf * H)

  beta_dt_day <- diag(wf) %*% (theta_day * xi[1]) %*% diag(1/as.vector(W_day))
  beta_dt_night <- diag(wf) %*% (theta_night * xi[2]) %*% diag(1/as.vector(W_night))
  beta_dt_manual <- (beta_dt_day * xi[1]) + (beta_dt_night * xi[2])

  # calculate fractional beta in Micro-MoB
  residency <- strata_to_residency_counts(H_counts = H_count)

  model <- setup_model_object()
  setup_human("strata", model = model, H = residency$H, J = residency$J)
  setup_timespent("dt", model = model, theta = theta)
  setup_biteweight("simple", model = model, wf = wf)

  Psi_t <- compute_Psi(human = model$human, xi = xi, t = 1)
  beta_dt <- compute_beta(human = model$human, Psi_t = Psi_t, t = 1)

  expect_equal(beta_dt[, , 1], beta_dt_day)
  expect_equal(beta_dt[, , 2], beta_dt_night)
  expect_equal(
    (beta_dt[, , 1] * xi[1]) + (beta_dt[, , 2] * xi[2]),
    beta_dt_manual
  )

  # test all beta matrices fulfill requirement on columns
  expect_true(all.equal(colSums(diag(H) %*% beta_dt_day), rep(1, 3)))
  expect_true(all.equal(colSums(diag(H) %*% beta_dt_night), rep(1, 3)))
  expect_true(all.equal(colSums(diag(H) %*% beta_dt_manual), rep(1, 3)))
  expect_true(all.equal(colSums(diag(H) %*% beta_dt[, , 1]), rep(1, 3)))
  expect_true(all.equal(colSums(diag(H) %*% beta_dt[, , 2]), rep(1, 3)))

  # bites are conserved
  bites <- pmax(1, rpois(n = p, lambda = 50))

  HBR_day <- beta_dt_day %*% bites
  bites_strata_day <- diag(H) %*% beta_dt_day %*% bites

  expect_equal(bites_strata_day / H, HBR_day)
  expect_equal((diag(H) %*% beta_dt[, , 1] %*% bites) / H, HBR_day)

  HBR_night <- beta_dt_night %*% bites
  bites_strata_night <- diag(H) %*% beta_dt_night %*% bites

  expect_equal(bites_strata_night / H, HBR_night)
  expect_equal((diag(H) %*% beta_dt[, , 2] %*% bites) / H, HBR_night)

  # daily computation is not expected to match when theta differs
  theta_average <- theta_day * xi[1] + theta_night * xi[2]

  model <- setup_model_object()
  setup_human("strata", model = model, H = residency$H, J = residency$J)
  setup_timespent("day", model = model, theta = theta_average)
  setup_biteweight("simple", model = model, wf = wf)

  Psi_t <- compute_Psi(human = model$human, xi = 1, t = 1)
  beta_day <- compute_beta(human = model$human, Psi_t = Psi_t, t = 1)

  expect_true(any(beta_dt_manual != beta_day))

  theta_average <- theta_day * 0.5 + theta_night * 0.5

  model <- setup_model_object()
  setup_human("strata", model = model, H = residency$H, J = residency$J)
  setup_timespent("day", model = model, theta = theta_average)
  setup_biteweight("simple", model = model, wf = wf)

  Psi_t <- compute_Psi(human = model$human, xi = 1, t = 1)
  beta_day <- compute_beta(human = model$human, Psi_t = Psi_t, t = 1)

  expect_true(any(beta_dt_manual != beta_day))
})


