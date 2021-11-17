test_that("test is_binary", {
  expect_true(is_binary(matrix(c(0, 1, 1, 0), 2, 2)))
  expect_true(is_binary(c(0, 1, 1, 0)))
  expect_false(is_binary(matrix(c(0, 12.5, 1, 0), 2, 2)))
  expect_false(is_binary(c(0, 1.25, 1, 0)))
})
