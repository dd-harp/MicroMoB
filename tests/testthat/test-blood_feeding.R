test_that("daily blood feeding test", {

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
  W <- t(theta) %*% (wf * as.vector(H))
  beta <- diag(wf) %*% theta %*% diag(1/as.vector(W))

  # calculate beta in Micro-MoB
  residency <- strata_to_residency_counts(H_counts = H_count)

  model <- new.env()
  setup_human("strata", model = model, H = residency$H, J = residency$J)
  setup_timespent("day", model = model, theta = theta)
  setup_biteweight("simple", model = model, wf = wf)

  beta_compare <- compute_beta(human = model$human, xi = 1, t = 1)

  expect_equal(beta, beta_compare)

})


test_that("fractional blood feeding test", {
  expect_equal(2 * 2, 4)
})
