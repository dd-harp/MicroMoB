test_that("strata to residency helper function", {

  J <- matrix(
    c(0.3, 0.5, 0.2,
      0.1, 0.6, 0.3), nrow = 3, ncol = 2, byrow = F
  )
  H <- c(50, 60)

  H_overall <- J %*% diag(H)

  residency <- strata_to_residency_proportion(H_strata = H, J_strata = J)

  expect_true(all(residency$H[residency$assignment_indices[1, ]] == H_overall[1, ]))
  expect_true(all(residency$H[residency$assignment_indices[2, ]] == H_overall[2, ]))
  expect_true(all(residency$H[residency$assignment_indices[3, ]] == H_overall[3, ]))

  expect_true(all(residency$H[residency$assignment_indices[, 1]] == H_overall[, 1]))
  expect_true(all(residency$H[residency$assignment_indices[, 2]] == H_overall[, 2]))

  expect_equal(residency$J %*% residency$H, J %*% H)

  # 2 strata, 3 patch
  J <- matrix(
    c(0.3, 0.5, 0.2,
      0.1, 0.6, 0.3), nrow = 3, ncol = 2, byrow = F
  )
  H <- c(s1 = 50, s2 = 60)

  bad_values <- c(Inf, NaN, NA, -5)
  for (v in bad_values) {
    J[1,2] <- v
    expect_error(strata_to_residency_proportion(H = H, J = J))
    expect_error(strata_to_residency_proportion(H = H, J = J))
    expect_error(strata_to_residency_proportion(H = H, J = J))
    expect_error(strata_to_residency_proportion(H = H, J = J))
  }

})


test_that("setting up human objects (strata) works when J is not specified", {

  expect_error(setup_human("strata", model = setup_model_object(), H = NULL))
  expect_error(setup_human("strata", model = setup_model_object(), H = rep(Inf, 5)))
  expect_error(setup_human("strata", model = setup_model_object(), H = rep(-5, 5)))
  expect_error(setup_human("strata", model = setup_model_object(), H = rep(NaN, 5)))
  expect_error(setup_human("strata", model = setup_model_object(), H = rep(NA, 5)))

  expect_error(setup_human("strata", model = setup_model_object(), H = c(Inf, 5)))
  expect_error(setup_human("strata", model = setup_model_object(), H = c(-5, 5)))
  expect_error(setup_human("strata", model = setup_model_object(), H = c(NaN, 5)))
  expect_error(setup_human("strata", model = setup_model_object(), H = c(NA, 5)))

  expect_error(setup_human("strata", model = setup_model_object(), H = c(1, 5), J = diag(5)))

  model = setup_model_object()
  setup_human("strata", model = model, H = c(1, 5))
  expect_equal(model$human$J, diag(2))
  expect_equal(model$human$H, c(1, 5))
  expect_true(inherits(model$human, "strata"))

  model = setup_model_object()
  expect_error(setup_human("blargh", model = model, H = c(1, 5)))
})


test_that("setting up human objects (strata) works when specifying J", {

  # 2 strata, 3 patch
  J_strata <- matrix(
    c(0.3, 0.5, 0.2,
      0.1, 0.6, 0.3), nrow = 3, ncol = 2, byrow = F
  )
  H_strata <- c(50, 60)

  residency <- strata_to_residency_proportion(H_strata = H_strata, J_strata = J_strata)

  J <- residency$J
  H <- residency$H
  n <- ncol(J)

  expect_error(setup_human("strata", model = setup_model_object(), H = NULL, J = J))
  expect_error(setup_human("strata", model = setup_model_object(), H = rep(Inf, n), J = J))
  expect_error(setup_human("strata", model = setup_model_object(), H = rep(-5, n), J = J))
  expect_error(setup_human("strata", model = setup_model_object(), H = rep(NaN, n), J = J))
  expect_error(setup_human("strata", model = setup_model_object(), H = rep(NA, n), J = J))
  expect_error(setup_human("strata", model = setup_model_object(), H = c(1, 2, 3), J = J))

  model = setup_model_object()
  setup_human("strata", model = model, H = H, J = J)

  expect_true(sum(model$human$J %*% model$human$H) == sum(H_strata))
  expect_equal(as.vector(J_strata %*% diag(H_strata)), model$human$H)

  human_pop <- matrix(
    c(15, 6,
      25, 36,
      10, 18), nrow = 3, ncol = 2, byrow = TRUE
  )
  expect_equal(as.vector(model$human$J %*% model$human$H), rowSums(human_pop))

  # test it works when we give J already formulated as binary
  model1 <- setup_model_object()
  setup_human("strata", model = model1, H = model$human$H, J = model$human$J)
  expect_true(all.equal(as.list(model$human), as.list(model1$human)))
})
