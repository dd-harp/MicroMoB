test_that("Test BQ deterministic dynamics (1-step)", {

  p <- 4
  l <- 5

  psiB <- c(0.25, 0.35, 0.4, 0.5)
  psiQ <- c(0.85, 0.75, 0.95, 0.90, 0.925)

  pB <- c(0.9, 0.95, 0.925, 0.91)
  pQ <- c(0.9, 0.95, 0.925, 0.91, 0.075)

  Psi_bb <- matrix(data = rexp(p*p), nrow = p, ncol = p)
  Psi_bb <- t(Psi_bb / rowSums(Psi_bb))
  Psi_bq <- matrix(data = rexp(p*l), nrow = p, ncol = l)
  Psi_bq <- t(Psi_bq / rowSums(Psi_bq))
  Psi_qb <- matrix(data = rexp(p*l), nrow = l, ncol = p)
  Psi_qb <- t(Psi_qb / rowSums(Psi_qb))
  Psi_qq <- matrix(data = rexp(l*l), nrow = l, ncol = l)
  Psi_qq <- t(Psi_qq / rowSums(Psi_qq))

  # kappa P(infection)
  kappa <- c(0.01, 0.05, 0.075, 0.1)

  # update matrices
  Mbq_inf <- Psi_bq %*% diag(pB*psiB*kappa, ncol = p)
  Mbq_noinf <- Psi_bq %*% diag(pB*psiB*(1-kappa), ncol = p)

  Mbb <- Psi_bb %*% diag(pB*(1-psiB), ncol = p)
  Mbq <- Psi_bq %*% diag(pB*psiB, ncol = p)
  Mqb <- Psi_qb %*% diag(pQ*psiQ, ncol = l)
  Mqq <- Psi_qq %*% diag(pQ*(1-psiQ), ncol = l)

  M <- rbind(
    cbind(Mbb, Mqb),
    cbind(Mbq, Mqq)
  )

  M_noinf <- rbind(
    cbind(Mbb, Mqb),
    cbind(Mbq_noinf, Mqq)
  )

  M_inf <- rbind(
    cbind(Mbb*0, Mqb*0),
    cbind(Mbq_inf, Mqq*0)
  )

  maxEIP <- 2

  EIP_shift <- matrix(data = 0, nrow = maxEIP + 1, ncol = maxEIP + 1)
  EIP_shift[1, 1] <- 1
  EIP_shift[2:(maxEIP+1), 1:maxEIP] <- diag(x = 1, nrow = maxEIP, ncol = maxEIP)

  M_mosy <- rpois(n = p+l, lambda = 1e4)
  Y_mosy <- matrix(data = rpois(n = (p+l)*(maxEIP+1), lambda = 1e3), nrow = p+l, ncol = maxEIP+1)

  q <- rep(1, p)
  f <- stats::qexp(p = psiB)

  # deterministic update
  mod <- make_MicroMoB(tmax = 2, p = p, l = l)
  setup_mosquito_BQ(model = mod, stochastic = FALSE, eip = maxEIP, pB = pB, pQ = pQ, psiQ = psiQ, Psi_bb = Psi_bb, Psi_bq = Psi_bq, Psi_qb = Psi_qb, Psi_qq = Psi_qq, M = M_mosy, Y = Y_mosy)
  setup_aqua_trace(model = mod, lambda = rep(0, l), stochastic = FALSE)

  newM <- M_noinf %*% M_mosy
  newY <- apply(X = Y_mosy, MARGIN = 2, FUN = function(Y) {M %*% Y})
  newY <- newY %*% EIP_shift
  newY[, 3] <- M_inf %*% M_mosy

  mod$mosquito$f <- f
  mod$mosquito$q <- q
  mod$mosquito$kappa <- kappa
  step_mosquitoes(model = mod)

  expect_equal(newM, mod$mosquito$M)
  expect_equal(newY, mod$mosquito$Y)

  # stochastic update
  mod <- make_MicroMoB(tmax = 2, p = p, l = l)
  setup_mosquito_BQ(model = mod, stochastic = TRUE, eip = maxEIP, pB = pB, pQ = pQ, psiQ = psiQ, Psi_bb = Psi_bb, Psi_bq = Psi_bq, Psi_qb = Psi_qb, Psi_qq = Psi_qq, M = M_mosy, Y = Y_mosy)
  setup_aqua_trace(model = mod, lambda = rep(0, l), stochastic = FALSE)

  mod$mosquito$f <- f
  mod$mosquito$q <- q
  mod$mosquito$kappa <- kappa
  step_mosquitoes(model = mod)

})
