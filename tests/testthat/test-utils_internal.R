test_that("test is_binary", {
  expect_true(is_binary(matrix(c(0, 1, 1, 0), 2, 2)))
  expect_true(is_binary(c(0, 1, 1, 0)))
  expect_false(is_binary(matrix(c(0, 12.5, 1, 0), 2, 2)))
  expect_false(is_binary(c(0, 1.25, 1, 0)))
})


test_that("distribute works", {
  p <- 2
  nn <- 5:10
  out <- vapply(X = nn, FUN = function(x) {distribute(n = x, p = p)}, FUN.VALUE = numeric(p))
  expect_equal(colSums(out), nn)

  p <- 3
  nn <- 5:10
  out <- vapply(X = nn, FUN = function(x) {distribute(n = x, p = p)}, FUN.VALUE = numeric(p))
  expect_equal(colSums(out), nn)
})


test_that("approx equal is working", {
  expect_true(all(approx_equal(1:5, 1:5)))
  expect_true(all(approx_equal(as.numeric(1:5), as.numeric(1:5))))
  expect_false(approx_equal(5, 5 + sqrt(.Machine$double.eps)))
})


test_that("time_patch_varying_parameter is working", {

  # test it errors with obviously bad input
  p <- 1
  tmax <- 5
  expect_error(time_patch_varying_parameter(param = numeric(0), p = p, tmax = tmax))
  expect_error(time_patch_varying_parameter(param = 1:20, p = p, tmax = tmax))
  expect_error(time_patch_varying_parameter(param = numeric(0), p = p, tmax = -1))
  expect_error(time_patch_varying_parameter(param = 1:20, p = p, tmax = NaN))
  expect_error(time_patch_varying_parameter(param = 1, p = -1, tmax = tmax))
  expect_error(time_patch_varying_parameter(param = 1, p = NaN, tmax = tmax))
  expect_error(time_patch_varying_parameter(param = rexp(90), p = p, tmax = tmax))

  p <- 3
  expect_error(time_patch_varying_parameter(param = matrix(data = rexp(35), nrow = 5, ncol = 7), p = p, tmax = tmax))
  expect_error(time_patch_varying_parameter(param = 1, p = p, tmax = -1))
  expect_error(time_patch_varying_parameter(param = 1, p = p, tmax = NaN))
  expect_error(time_patch_varying_parameter(param = 1, p = -1, tmax = tmax))
  expect_error(time_patch_varying_parameter(param = 1, p = NaN, tmax = tmax))

  # works with p = 1
  p <- 1
  tmax <- 10
  const_vec <- rexp(n = p)
  period_vec <- rexp(n = 365)
  full_vec <- rexp(n = tmax)
  period_mat <- matrix(rexp(n = 365), nrow = p, ncol = 365)
  full_mat <- matrix(rexp(n = p*tmax), nrow = p, ncol = tmax)

  const_vec_out <- time_patch_varying_parameter(param = const_vec, p = p, tmax = tmax)
  expect_equal(nrow(const_vec_out), p)
  expect_equal(ncol(const_vec_out), tmax)
  expect_equal(as.vector(const_vec_out), rep(const_vec, tmax))

  period_vec_out <- time_patch_varying_parameter(param = period_vec, p = p, tmax = tmax)
  expect_equal(nrow(period_vec_out), p)
  expect_equal(ncol(period_vec_out), tmax)
  expect_equal(as.vector(period_vec_out), period_vec[1:tmax])

  full_vec_out <- time_patch_varying_parameter(param = full_vec, p = p, tmax = tmax)
  expect_equal(nrow(full_vec_out), p)
  expect_equal(ncol(full_vec_out), tmax)
  expect_equal(as.vector(full_vec_out), full_vec)

  period_mat_out <- time_patch_varying_parameter(param = period_mat, p = p, tmax = tmax)
  expect_equal(nrow(period_mat_out), p)
  expect_equal(ncol(period_mat_out), tmax)
  expect_equal(period_mat_out, period_mat[, 1:tmax, drop = FALSE])

  full_mat_out <- time_patch_varying_parameter(param = full_mat, p = p, tmax = tmax)
  expect_equal(nrow(full_mat_out), p)
  expect_equal(ncol(full_mat_out), tmax)
  expect_equal(full_mat_out, full_mat)

  # works with p > 1
  p <- 3
  tmax <- 10
  const_vec <- rexp(n = p)
  period_vec <- rexp(n = 365)
  full_vec <- rexp(n = tmax)
  period_mat <- matrix(rexp(n = 365), nrow = p, ncol = 365)
  full_mat <- matrix(rexp(n = p*tmax), nrow = p, ncol = tmax)

  const_vec_out <- time_patch_varying_parameter(param = const_vec, p = p, tmax = tmax)
  expect_equal(nrow(const_vec_out), p)
  expect_equal(ncol(const_vec_out), tmax)
  expect_equal(const_vec_out, replicate(n = tmax, expr = const_vec))

  period_vec_out <- time_patch_varying_parameter(param = period_vec, p = p, tmax = tmax)
  expect_equal(nrow(period_vec_out), p)
  expect_equal(ncol(period_vec_out), tmax)
  expect_equal(period_vec_out, t(replicate(n = p, expr = period_vec[1:tmax])))

  full_vec_out <- time_patch_varying_parameter(param = full_vec, p = p, tmax = tmax)
  expect_equal(nrow(full_vec_out), p)
  expect_equal(ncol(full_vec_out), tmax)
  expect_equal(full_vec_out, t(replicate(n = p, expr = full_vec)))

  period_mat_out <- time_patch_varying_parameter(param = period_mat, p = p, tmax = tmax)
  expect_equal(nrow(period_mat_out), p)
  expect_equal(ncol(period_mat_out), tmax)
  expect_equal(period_mat_out, period_mat[, 1:tmax, drop = FALSE])

  full_mat_out <- time_patch_varying_parameter(param = full_mat, p = p, tmax = tmax)
  expect_equal(nrow(full_mat_out), p)
  expect_equal(ncol(full_mat_out), tmax)
  expect_equal(full_mat_out, full_mat)

})


test_that("time_patch_varying_parameter is working", {

  # test it errors with obviously bad input
  tmax <- 5
  expect_error(time_varying_parameter(param = numeric(0), tmax = tmax))
  expect_error(time_varying_parameter(param = 1:20, tmax = tmax))
  expect_error(time_varying_parameter(param = numeric(0), tmax = -1))
  expect_error(time_varying_parameter(param = 1:20, tmax = NaN))
  expect_error(time_varying_parameter(param = rexp(90), tmax = tmax))

  # works with p = 1
  tmax <- 10
  const <- rexp(n = 1)
  period <- rexp(n = 365)
  full <- rexp(n = tmax)

  const_out <- time_varying_parameter(param = const, tmax = tmax)
  expect_equal(length(const_out), tmax)
  expect_equal(const_out, rep(const, tmax))

  period_out <- time_varying_parameter(param = period, tmax = tmax)
  expect_equal(length(period_out), tmax)
  expect_equal(period_out, period[1:tmax])

  full_out <- time_varying_parameter(param = full, tmax = tmax)
  expect_equal(length(full_out), tmax)
  expect_equal(full_out, full)

})

