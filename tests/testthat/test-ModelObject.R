test_that("the model object constructs", {
  obj <- make_MicroMoB(tmax = 10, p = 1)
  expect_true(inherits(obj, "MicroMoB"))
  expect_true(inherits(obj$global, "list"))
  expect_true(storage.mode(obj) == "environment")
  expect_equal(obj$global$tmax, 10)
  expect_equal(obj$global$tnow, 1)
  expect_equal(obj$global$p, obj$global$l)

  expect_error(make_MicroMoB(tmax = Inf, p = 1))
  expect_error(make_MicroMoB(tmax = -5, p = 1))
  expect_error(make_MicroMoB(tmax = 0, p = 1))
})
