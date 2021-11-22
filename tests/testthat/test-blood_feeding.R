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

  wf <- c(rep(0.9, p), rep(1.1, p))

  # calculate beta manually
  W <- t(theta) %*% (wf * H)
  beta <- diag(wf) %*% theta %*% diag(1/as.vector(W))
  beta <- beta %*% diag(1/colSums(beta))

  # calculate beta in Micro-MoB
  residency <- strata_to_residency_counts(H_counts = H_count)

  model <- new.env()
  setup_human("strata", model = model, H = residency$H, J = residency$J)
  setup_timespent("day", model = model, theta = theta)
  setup_biteweight("simple", model = model, wf = wf)

  beta_compare <- compute_beta(human = model$human, xi = 1, t = 1)

  expect_equal(beta, beta_compare)

})


test_that("fractional beta, biting distribution matrix test", {

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
  theta_night[, 1] <- theta_night[, 1] * 0.8
  theta_night[, 2] <- theta_night[, 2] * 0.5
  theta_night[, 3] <- theta_night[, 3] * 1.5
  theta_night <- theta_night / rowSums(theta_night)

  theta <- list(theta_day, theta_night)

  wf <- c(rep(0.9, p), rep(1.1, p))
  xi <- c(0.25, 0.75)

  # calculate beta manually
  W_day <- (t(theta_day) * xi[1]) %*% (wf * H)
  W_night <- (t(theta_night) * xi[2]) %*% (wf * H)

  beta_day <- diag(wf) %*% (theta_day * xi[1]) %*% diag(1/as.vector(W_day))
  beta_day <- beta_day %*% diag(1/colSums(beta_day))
  beta_night <- diag(wf) %*% (theta_night * xi[2]) %*% diag(1/as.vector(W_night))
  beta_night <- beta_night %*% diag(1/colSums(beta_night))

  # calculate beta in Micro-MoB
  residency <- strata_to_residency_counts(H_counts = H_count)

  model <- new.env()
  setup_human("strata", model = model, H = residency$H, J = residency$J)
  setup_timespent("dt", model = model, theta = theta)
  setup_biteweight("simple", model = model, wf = wf)

  beta_compare <- compute_beta(human = model$human, xi = xi, t = 1)

  expect_equal(beta_day, beta_compare[, , 1])
  expect_equal(beta_night, beta_compare[, , 2])

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
  beta_dt_day <- beta_dt_day %*% diag(1/colSums(beta_dt_day))
  beta_dt_night <- diag(wf) %*% (theta_night * xi[2]) %*% diag(1/as.vector(W_night))
  beta_dt_night <- beta_dt_night %*% diag(1/colSums(beta_dt_night))
  beta_dt_manual <- (beta_dt_day * xi[1]) + (beta_dt_night * xi[2])

  # calculate fractional beta in Micro-MoB
  residency <- strata_to_residency_counts(H_counts = H_count)

  model <- new.env()
  setup_human("strata", model = model, H = residency$H, J = residency$J)
  setup_timespent("dt", model = model, theta = theta)
  setup_biteweight("simple", model = model, wf = wf)

  beta_dt <- compute_beta(human = model$human, xi = xi, t = 1)

  expect_equal(beta_dt[, , 1], beta_dt_day)
  expect_equal(beta_dt[, , 2], beta_dt_night)

  # calculate daily beta
  model <- new.env()
  setup_human("strata", model = model, H = residency$H, J = residency$J)
  setup_timespent("day", model = model, theta = theta[[1]])
  setup_biteweight("simple", model = model, wf = wf)

  beta_daily <- compute_beta(human = model$human, t = 1)

  # calculate fractional  beta manually
  theta_daily <- theta_day + theta_night
  theta_daily <- diag(1/rowSums(theta_daily)) %*% theta_daily

  W <- t(theta_daily) %*% (wf * H)
  beta_daily_manual <- diag(wf) %*% theta_daily %*% diag(1/as.vector(W))
  beta_daily_manual <- beta_daily_manual %*% diag(1/colSums(beta_daily_manual))

  expect_equal(beta_dt_manual, beta_daily)
  expect_equal(beta_daily_manual, beta_daily)

})



test_that("daily v. fractional (3 chunks) beta, biting distribution matrix test", {
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


})


