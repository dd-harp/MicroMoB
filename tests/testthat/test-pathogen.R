test_that("RM pathogen works", {
  model <- new.env()
  setup_human("strata", model = model, H = c(10, 50, 20))
  expect_error(setup_pathogen("blargh", model = model))
  expect_error(setup_pathogen("rm", model = model))
  expect_error(setup_pathogen("rm", model = model, rm_parameters = list(zzz = 5, X = c(0.5, 0.2, 0.5))))

  rm_parameters <- list(X = c(0.5, 0.2, 0.5), r = 1/200, b = 0.55, c = 0.15)
  setup_pathogen("rm", model = model, rm_parameters = rm_parameters)

  expect_true(all(compute_x(pathogen = model$pathogen) == rm_parameters$X * rm_parameters$c))
  expect_equal(model$pathogen$r, rm_parameters$r)
  expect_equal(model$pathogen$b, rm_parameters$b)
  expect_equal(model$pathogen$c, rm_parameters$c)

})
