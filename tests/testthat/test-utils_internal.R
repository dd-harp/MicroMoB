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


test_that("time varying parameter working", {
  p <- 1
  tmax <- 10

  p <- 1
  tmax <- 365

  p <- 1
  tmax <- 400
})
