test_that("strata to residency helper function", {

  J <- matrix(
    c(0.3, 0.5, 0.2,
      0.1, 0.6, 0.3), nrow = 3, ncol = 2, byrow = F
  )
  H <- c(50, 60)

  H_overall <- J %*% diag(H)

  residency <- strata_to_residency(H_strata = H, J_strata = J)

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
    expect_error(strata_to_residency(H = H, J = J))
    expect_error(strata_to_residency(H = H, J = J))
    expect_error(strata_to_residency(H = H, J = J))
    expect_error(strata_to_residency(H = H, J = J))
  }

})


test_that("setting up human objects (strata) works when J is not specified", {

  expect_error(setup.human("strata", model = new.env(), H = NULL))
  expect_error(setup.human("strata", model = new.env(), H = rep(Inf, 5)))
  expect_error(setup.human("strata", model = new.env(), H = rep(-5, 5)))
  expect_error(setup.human("strata", model = new.env(), H = rep(NaN, 5)))
  expect_error(setup.human("strata", model = new.env(), H = rep(NA, 5)))

  expect_error(setup.human("strata", model = new.env(), H = c(Inf, 5)))
  expect_error(setup.human("strata", model = new.env(), H = c(-5, 5)))
  expect_error(setup.human("strata", model = new.env(), H = c(NaN, 5)))
  expect_error(setup.human("strata", model = new.env(), H = c(NA, 5)))

  expect_error(setup.human("strata", model = new.env(), H = c(1, 5), J = diag(5)))

  model <- new.env()
  setup.human("strata", model = model, H = c(1, 5))
  expect_equal(model$human$J, diag(2))
  expect_equal(model$human$H, c(1, 5))
  expect_true(inherits(model$human, "strata"))
})


test_that("setting up human objects (strata) works when specifying J", {

  # 2 strata, 3 patch
  J_strata <- matrix(
    c(0.3, 0.5, 0.2,
      0.1, 0.6, 0.3), nrow = 3, ncol = 2, byrow = F
  )
  H_strata <- c(50, 60)

  residency <- strata_to_residency(H_strata = H_strata, J_strata = J_strata)

  J <- residency$J
  H <- residency$H
  n <- ncol(J)

  expect_error(setup.human("strata", model = new.env(), H = NULL, J = J))
  expect_error(setup.human("strata", model = new.env(), H = rep(Inf, n), J = J))
  expect_error(setup.human("strata", model = new.env(), H = rep(-5, n), J = J))
  expect_error(setup.human("strata", model = new.env(), H = rep(NaN, n), J = J))
  expect_error(setup.human("strata", model = new.env(), H = rep(NA, n), J = J))
  expect_error(setup.human("strata", model = new.env(), H = c(1, 2, 3), J = J))

  model <- new.env()
  setup.human("strata", model = model, H = H, J = J)

  expect_true(sum(model$human$J %*% model$human$H) == sum(H_strata))
  expect_equal(as.vector(J_strata %*% diag(H_strata)), model$human$H)

  human_pop <- matrix(
    c(15, 6,
      25, 36,
      10, 18), nrow = 3, ncol = 2, byrow = TRUE
  )
  expect_equal(as.vector(model$human$J %*% model$human$H), rowSums(human_pop))

  # test it works when we give J already formulated as binary
  model1 <- new.env()
  setup.human("strata", model = model1, H = model$human$H, J = model$human$J)
  expect_true(all.equal(as.list(model), as.list(model1)))
})


test_that("setting up time spent (day) works", {

  model <- new.env()
  setup.human("strata", model = model, H = rep(10, 3))
  setup.timespent("day", model = model)

  expect_true(inherits(model$tisp, "day"))
  expect_equal(model$tisp$Theta_t, diag(3))

  # 2 patch, 3 strata
  H <- c(50, 20, 10, 30, 20, 10)
  J <- matrix(
    c(1, 0, 1, 0, 1, 0,
      0, 1, 0, 1, 0, 1),
    nrow = 2, ncol = 6, byrow = TRUE
  )

  # theta keeps everyone in one place
  model <- new.env()
  setup.human("strata", model = model, H = H, J = J)
  setup.timespent("day", model = model)

  expect_equal(model$tisp$Theta_t %*% H, model$human$J %*% model$human$H)
  expect_equal(sum(model$tisp$Theta_t %*% H), sum(H))
})

test_that("setting up time spent (dt) works", {

  # 2 patch, 3 strata
  H <- c(50, 20, 10, 30, 20, 10)
  J <- matrix(
    c(1, 0, 1, 0, 1, 0,
      0, 1, 0, 1, 0, 1),
    nrow = 2, ncol = 6, byrow = TRUE
  )

  expect_equal(as.vector(J %*% H), c(sum(H[1], H[3], H[5]), sum(H[2], H[4], H[6])))

  theta_day <- matrix(
    c(0.4, 0.6,
      0.5, 0.5,
      0.2, 0.8,
      0.5, 0.5,
      0.6, 0.4,
      0.2, 0.8), nrow = length(H), ncol = 2, byrow = TRUE
  )

  theta_night <- matrix(
    c(0.9, 0.1,
      0.6, 0.4,
      0.5, 0.5,
      0.6, 0.4,
      0.8, 0.2,
      0.4, 0.6), nrow = length(H), ncol = 2, byrow = TRUE
  )

  theta <- list(theta_day, theta_night)

  model <- new.env()
  setup.human("strata", model = model, H = H, J = J)
  setup.timespent("dt", model = model, theta = theta)
  setup.biteweight("simple", model = model)

  xi <- c(0.15, 0.85)
  model$tisp$Psi_t <- compute_Psi.timespent(tisp = model$tisp, xi = xi, t = 1)
  W <- compute_W.timespent(tisp = model$tisp, biteweight = model$biteweight, human = model$human, t = 1)

  wt <- rep(1, length(H))
  W_expected <- ((t(theta_day) %*% (wt * H)) * (xi[1])) + ((t(theta_night) %*% (wt * H)) * (xi[2]))

  expect_equal(W, W_expected)

  # 3 strata, 2 patch
  H_strata <- c(70, 40, 30)
  J_strata <- matrix(
    c(50/70,  10/40, 20/30,
      1-50/70, 1-10/40, 1-20/30),
    nrow = 2,
    ncol = 3, byrow = TRUE
  )

  residency <- strata_to_residency(H_strata = H_strata, J_strata = J_strata)
  J <- residency$J
  H <- residency$H

  model <- new.env()
  setup.human("strata", model = model, H = H, J = J)
  setup.timespent("dt", model = model, theta = theta)
  setup.biteweight("simple", model = model)

  xi <- c(0.15, 0.85)
  model$tisp$Psi_t <- compute_Psi.timespent(tisp = model$tisp, xi = xi, t = 1)
  W <- compute_W.timespent(tisp = model$tisp, biteweight = model$biteweight, human = model$human, t = 1)

  expect_equal(W, W_expected)
})
