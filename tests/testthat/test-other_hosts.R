test_that("simple other hosts works correctly", {
  model <- setup_model_object()
  setup_human("strata", model = model, H = rep(10, 3))
  expect_error(setup_otherhosts("blargh", model))

  setup_otherhosts("simple", model)
  expect_true(!is.null(model$otherhosts))
  expect_true(inherits(model$otherhosts, "simple"))
  expect_equal(compute_B(model$otherhosts), rep(0, 3))

  class(model$otherhosts) <- "blargh"
  expect_error(compute_B(model$otherhosts))
})


test_that("simple visitors works correctly", {
  model <- setup_model_object()
  setup_human("strata", model = model, H = rep(10, 3))
  expect_error(setup_visitors("blargh", model))

  setup_visitors("simple", model)
  expect_true(!is.null(model$visitors))
  expect_true(inherits(model$visitors, "simple"))
  expect_equal(compute_W_delta(model$visitors), rep(0, 3))
  expect_equal(compute_x_delta(model$visitors), rep(0, 3))

  class(model$visitors) <- "blargh"
  expect_error(compute_B(model$visitors))

  setup_visitors("simple", model, W_delta = c(0.5, 0.25, 0.3), x_delta = c(0.7, 0.2, 0.1))
  expect_equal(compute_W_delta(model$visitors), c(0.5, 0.25, 0.3))
  expect_equal(compute_x_delta(model$visitors), c(0.7, 0.2, 0.1))
})
