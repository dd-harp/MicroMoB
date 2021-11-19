test_that("simple other hosts works correctly", {
  model <- new.env()
  setup.human("strata", model = model, H = rep(10, 3))
  expect_error(setup_otherhosts("blargh", model))

  setup_otherhosts("simple", model)
  expect_true(!is.null(model$otherhosts))
  expect_true(inherits(model$otherhosts, "simple"))
  expect_equal(compute_B(model$otherhosts), rep(0, 3))

  class(model$otherhosts) <- "blargh"
  expect_error(compute_B(model$otherhosts))
})
