test_that("sample draw_multinomial", {

  n <- 10
  prob <- c(0.5, 0.2, 0.05, 0.25)
  draw <- draw_multinom(n = n, prob = prob)
  expect_equal(sum(draw), n)
  expect_length(draw, length(prob))

})

test_that("sample stochastic vector", {

  cell <- 4
  bins <- 4
  pmat <- matrix(rexp(cell*bins), nrow = cell, ncol = bins)
  pmat <- pmat / rowSums(pmat)

  x <- rpois(n = cell, lambda = 1e5)

  expect <- as.vector(x %*% pmat)

  draw <- sample_stochastic_vector(x = x, prob = pmat)
  expect_true(all(abs(expect - draw) / expect < 0.05))
  expect_equal(sum(x), sum(draw))

  cell <- 4
  bins <- 5
  pmat <- matrix(rexp(cell*bins), nrow = cell, ncol = bins)
  pmat <- pmat / rowSums(pmat)

  x <- rpois(n = cell, lambda = 1e5)

  expect <- as.vector(x %*% pmat)

  draw <- sample_stochastic_vector(x = x, prob = pmat)
  expect_true(all(abs(expect - draw) / expect < 0.05))
  expect_equal(sum(x), sum(draw))

})


test_that("sample stochastic matrix", {
  ZZ <- matrix(data = rpois(n = 2*3, lambda = (1:6) * 1e7),nrow = 2,ncol = 3)
  psi <- matrix(c(
    0.9, 0.025, 0.075,
    0.1, 0.7, 0.2,
    0.01, 0.09, 0.9
  ),nrow=3,ncol=3,byrow=T)

  expectation <- ZZ %*% psi

  draw <- sample_stochastic_matrix(x = ZZ, prob = psi)

  expect_true(all(abs(expectation - draw) / expectation < 0.05))
  expect_equal(sum(draw), sum(ZZ))
  expect_equal(rowSums(draw), rowSums(ZZ))
})


test_that("proportional strata to residency functions working", {

  # p > s case
  s <- 2
  p <- 3
  n <- s*p

  J <- matrix(
    data = c(
      0.5, 0.6,
      0.2, 0.3,
      0.3, 0.1
    ), nrow = p, ncol = s, byrow = TRUE
  )

  H <- c(100, 120)

  H_expect <- c(50, 20, 30, 72, 36, 12)
  J_expect <- cbind(diag(3), diag(3))

  J_out <- strata_to_residency_proportion(H_strata = H, J_strata = J)

  expect_equal(nrow(J_out$J), p)
  expect_equal(ncol(J_out$J), n)
  expect_equal(J_out$J, J_expect)
  expect_equal(J_out$H, H_expect)
  expect_equal(as.vector(J %*% diag(H)), J_out$H)

  # s > p case
  s <- 3
  p <- 2
  n <- s*p

  J <- matrix(
    data = c(
      0.4, 0.2, 0.9,
      0.6, 0.8, 0.1
    ), nrow = p, ncol = s, byrow = TRUE
  )

  H <- c(100, 120, 140)

  H_expect <- c(40, 60, 24, 96, 126, 14)
  J_expect <- do.call(cbind, replicate(n = 3, expr = diag(2),simplify = FALSE))

  J_out <- strata_to_residency_proportion(H_strata = H, J_strata = J)

  expect_equal(nrow(J_out$J), p)
  expect_equal(ncol(J_out$J), n)
  expect_equal(J_out$J, J_expect)
  expect_equal(J_out$H, H_expect)
  expect_equal(as.vector(J %*% diag(H)), J_out$H)

})


test_that("counts strata to residency functions working", {

  # p > s case
  s <- 2
  p <- 3
  n <- s*p

  H_counts <- matrix(
    data = c(
      20, 30,
      50, 20,
      10, 35
    ), nrow = p, ncol = s, byrow = TRUE
  )

  J_expect <- cbind(diag(3), diag(3))
  H_expect <- as.vector(H_counts)

  J_out <- strata_to_residency_counts(H_counts = H_counts)

  expect_equal(nrow(J_out$J), p)
  expect_equal(ncol(J_out$J), n)
  expect_equal(J_out$J, J_expect)
  expect_equal(J_out$H, H_expect)
  expect_equal(as.vector(J_out$J %*% J_out$H), rowSums(H_counts))

  # s > p case
  s <- 3
  p <- 2
  n <- s*p

  H_counts <- matrix(
    data = c(
      20, 60, 10,
      30, 85, 25
    ), nrow = p, ncol = s, byrow = TRUE
  )

  J_expect <- cbind(diag(2), diag(2), diag(2))
  H_expect <- as.vector(H_counts)

  J_out <- strata_to_residency_counts(H_counts = H_counts)

  expect_equal(nrow(J_out$J), p)
  expect_equal(ncol(J_out$J), n)
  expect_equal(J_out$J, J_expect)
  expect_equal(J_out$H, H_expect)
  expect_equal(as.vector(J_out$J %*% J_out$H), rowSums(H_counts))

})
