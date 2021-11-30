test_that("RM pathogen works", {
  model <- setup_model_object()
  H <- c(10, 50, 20)
  setup_human("strata", model = model, H = H)
  expect_error(setup_pathogen("blargh", model = model))
  expect_error(setup_pathogen("rm", model = model))
  expect_error(setup_pathogen("rm", model = model, rm_parameters = list(zzz = 5, X = c(0.5, 0.2, 0.5), stochastic = FALSE)))

  rm_parameters <- list(X = c(0.5, 0.2, 0.5), r = 1/200, b = 0.55, c = 0.15, stochastic = FALSE)
  setup_pathogen("rm", model = model, rm_parameters = rm_parameters)

  expect_true(all(compute_x(model = model) == (rm_parameters$X / H) * rm_parameters$c))
  expect_equal(model$pathogen$r, rep(rm_parameters$r, length(H)))
  expect_equal(model$pathogen$b, diag(rep(rm_parameters$b, length(H))))
  expect_equal(model$pathogen$c, rep(rm_parameters$c, length(H)))

})
