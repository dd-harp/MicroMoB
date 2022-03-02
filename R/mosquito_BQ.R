# classes and methods to implement a simple behavioral state model of mosquitoes with infection and incubation

#' @title Setup blood feeding & oviposition (BQ) behavioral state mosquito model
#' @description This is a behavioral state model which allows for time varying EIP and
#' survival probability. Mosquitoes transition between blood feeding (B) and
#' oviposition (Q) depending on teh success (or not) of those biological activities.
#' It complies with the mosquito component interface, and may be simulated deterministically or stochastically.
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @param stochastic should the model update deterministically or stochastically?
#' @param eip the Extrinsic Incubation Period (may be time varying see [MicroMoB::time_varying_parameter])
#' @param pB daily survival probability during blood feeding (may be time and patch varying see [MicroMoB::time_patch_varying_parameter])
#' @param pQ daily survival probability during oviposition (may be time and patch varying see [MicroMoB::time_patch_varying_parameter])
#' @param psiQ oviposition success probability (may be time and patch varying see [MicroMoB::time_patch_varying_parameter])
#' @param Psi_bb movement matrix from blood feeding haunts to blood feeding haunts (columns must sum to 1, `p` rows and columns)
#' @param Psi_bq movement matrix from blood feeding haunts to aquatic habitats (columns must sum to 1, `l` rows and `p` columns)
#' @param Psi_qb movement matrix from aquatic habitats to blood feeding haunts (columns must sum to 1, `p` rows and `l` columns)
#' @param Psi_qq movement matrix from aquatic habitats to aquatic habitats (columns must sum to 1, `l` rows and columns)
#' @param nu number of eggs laid per oviposition
#' @param M number of susceptible mosquitoes (vector of length `p + l`)
#' @param Y number of incubating mosquitoes (vector of length `p + l`)
#' @param Z number of infectious mosquitoes (vector of length `p + l`)
#' @return no return value
#' @export
setup_mosquito_BQ <- function(model, stochastic, eip, pB, pQ, psiQ, Psi_bb, Psi_bq, Psi_qb, Psi_qq, nu = 25, M, Y) {
  stopifnot(inherits(model, "MicroMoB"))
  stopifnot(is.logical(stochastic))

  p <- model$global$p
  l <- model$global$l
  tmax <- model$global$tmax

  eip_vec <- time_varying_parameter(param = eip, tmax = tmax)

  maxEIP <- max(eip_vec)

  pB_mat <- time_patch_varying_parameter(param = pB, p = p, tmax = tmax)
  pQ_mat <- time_patch_varying_parameter(param = pQ, p = l, tmax = tmax)
  psiQ_mat <- time_patch_varying_parameter(param = psiQ, p = l, tmax = tmax)

  stopifnot(dim(Psi_bb) == p)
  stopifnot(dim(Psi_bq) == c(l, p))
  stopifnot(dim(Psi_qb) == c(p, l))
  stopifnot(dim(Psi_qq) == l)
  stopifnot(approx_equal(colSums(Psi_bb), 1))
  stopifnot(approx_equal(colSums(Psi_bq), 1))
  stopifnot(approx_equal(colSums(Psi_qb), 1))
  stopifnot(approx_equal(colSums(Psi_qq), 1))

  stopifnot(length(nu) == 1L)

  mosy_class <- c("BQ")
  if (stochastic) {
    mosy_class <- c(mosy_class, "BQ_stochastic")
  } else {
    mosy_class <- c(mosy_class, "BQ_deterministic")
  }

  model$mosquito <- structure(list(), class = mosy_class)
  model$mosquito$psiQ_mat <- psiQ_mat
  model$mosquito$nu <- nu
  model$mosquito$eip <- eip_vec
  model$mosquito$maxEIP <- maxEIP

  model$mosquito$pB_mat <- pB_mat
  model$mosquito$pQ_mat <- pQ_mat

  model$mosquito$Psi_bb <- Psi_bb
  model$mosquito$Psi_bq <- Psi_bq
  model$mosquito$Psi_qb <- Psi_qb
  model$mosquito$Psi_qq <- Psi_qq

  model$mosquito$kappa <- rep(0, p)
  model$mosquito$f <- rep(0, p)
  model$mosquito$q <- rep(0, p)

  stopifnot(length(M) == p+l)
  stopifnot(nrow(Y) == p+l)
  stopifnot(ncol(Y) == maxEIP+1)
  stopifnot(is.finite(M))
  stopifnot(is.finite(Y))
  stopifnot(M >= 0)
  stopifnot(Y >= 0)

  if (stochastic) {
    model$mosquito$M <- matrix(data = M, ncol = 1)
  } else {
    model$mosquito$M <- M
  }
  model$mosquito$Y <- Y

  # matrix which multiplies Y on the right to shift all by one day
  EIP_shift <- matrix(data = 0, nrow = maxEIP + 1, ncol = maxEIP + 1)
  EIP_shift[1, 1] <- 1
  EIP_shift[2:(maxEIP+1), 1:maxEIP] <- diag(x = 1, nrow = maxEIP, ncol = maxEIP)
  model$mosquito$EIP_shift <- EIP_shift

}


#


# update mosquitoes over one time step

#' @title Update blood feeding & oviposition (BQ) behavioral state mosquitoes
#' @description This function dispatches on the second argument of `model$mosquito`
#' for stochastic or deterministic behavior.
#' @inheritParams step_mosquitoes
#' @return no return value
#' @details see [MicroMoB::step_mosquitoes.BQ_deterministic] and [MicroMoB::step_mosquitoes.BQ_stochastic]
#' @export
step_mosquitoes.BQ <- function(model) {
  NextMethod()
}

#' @title Update blood feeding & oviposition (BQ) behavioral state mosquitoes (deterministic)
#' @inheritParams step_mosquitoes
#' @return no return value
#' @importFrom stats pexp
#' @export
step_mosquitoes.BQ_deterministic <- function(model) {

  # parameters
  tnow <- model$global$tnow
  EIP <- model$mosquito$eip[tnow]
  maxEIP <- model$mosquito$maxEIP
  p <- model$global$p
  l <- model$global$l

  psiB <- pexp(q = model$mosquito$f * model$mosquito$q)
  psiQ <- model$mosquito$psiQ_mat[, tnow]

  pB <- model$mosquito$pB_mat[, tnow]
  pQ <- model$mosquito$pQ_mat[, tnow]

  kappa <- model$mosquito$kappa

  # construct update matrices
  Mbq_inf <- model$mosquito$Psi_bq %*% diag(x = pB*psiB*kappa, ncol = p)
  Mbq_noinf <- model$mosquito$Psi_bq %*% diag(x = pB*psiB*(1-kappa), ncol = p)

  Mbb <- model$mosquito$Psi_bb %*% diag(x = pB*(1-psiB), ncol = p)
  Mbq <- model$mosquito$Psi_bq %*% diag(x = pB*psiB, ncol = p)
  Mqb <- model$mosquito$Psi_qb %*% diag(x = pQ*psiQ, ncol = l)
  Mqq <- model$mosquito$Psi_qq %*% diag(x = pQ*(1-psiQ), ncol = l)

  M <- matrix(data = 0, nrow = l+p, ncol = l+p)
  M[1:p, 1:p] <- Mbb
  M[1:p, (p+1):(l+p)] <- Mqb
  M[(p+1):(l+p), 1:p] <- Mbq
  M[(p+1):(l+p), (p+1):(l+p)] <- Mqq

  M_noinf <- matrix(data = 0, nrow = l+p, ncol = l+p)
  M_noinf[1:p, 1:p] <- Mbb
  M_noinf[1:p, (p+1):(l+p)] <- Mqb
  M_noinf[(p+1):(l+p), 1:p] <- Mbq_noinf
  M_noinf[(p+1):(l+p), (p+1):(l+p)] <- Mqq

  M_inf <- matrix(data = 0, nrow = l+p, ncol = l+p)
  M_inf[(p+1):(l+p), 1:p] <- Mbq_inf

  # update
  Y0 <- M_inf %*% model$mosquito$M
  model$mosquito$M <- M_noinf %*% model$mosquito$M
  for (i in 1:(maxEIP+1)) {
    model$mosquito$Y[, i] <- M %*% model$mosquito$Y[, i]
  }
  model$mosquito$Y <- model$mosquito$Y %*% model$mosquito$EIP_shift
  model$mosquito$Y[, EIP+1L] <- Y0

  # newly emerging adults
  lambda <- compute_emergents(model)
  model$mosquito$M[1:p, 1] <- model$mosquito$M[1:p, 1] + model$mosquito$Psi_qb %*% lambda

}


#' @title Update blood feeding & oviposition (BQ) behavioral state mosquitoes (stochastic)
#' @inheritParams step_mosquitoes
#' @return no return value
#' @importFrom stats pexp rbinom
#' @export
step_mosquitoes.BQ_stochastic <- function(model) {

  # parameters
  tnow <- model$global$tnow
  EIP <- model$mosquito$eip[tnow]
  maxEIP <- model$mosquito$maxEIP
  p <- model$global$p
  l <- model$global$l

  psiB <- pexp(q = model$mosquito$f * model$mosquito$q)
  psiQ <- model$mosquito$psiQ_mat[, tnow]

  pB <- model$mosquito$pB_mat[, tnow]
  pQ <- model$mosquito$pQ_mat[, tnow]

  kappa <- model$mosquito$kappa

  # create movement matrices for success and fail
  Psi_success <- matrix(data = 0, nrow = l+p, ncol = l+p)
  Psi_success[1:p, (p+1):(l+p)] <- model$mosquito$Psi_qb
  Psi_success[(p+1):(l+p), 1:p] <- model$mosquito$Psi_bq

  Psi_fail <- matrix(data = 0, nrow = l+p, ncol = l+p)
  Psi_fail[1:p, 1:p] <- model$mosquito$Psi_bb
  Psi_fail[(p+1):(l+p), (p+1):(l+p)] <- model$mosquito$Psi_qq

  # sample survivors
  model$mosquito$M <- rbinom(n = l+p, size = model$mosquito$M, prob = c(pB, pQ))
  for (i in 1:(maxEIP+1)) {
    model$mosquito$Y[, i] <- rbinom(n = l+p, size = model$mosquito$Y[, i], prob = c(pB, pQ))
  }

  # sample behavioral success and movement
  for (i in 1:(maxEIP+1)) {
    BQ_success <- rbinom(n = l+p, size = model$mosquito$Y[, i], prob = c(psiB, psiQ))
    BQ_move <- sample_stochastic_vector(x = BQ_success, prob = t(Psi_success))
    BQ_move <- BQ_move + sample_stochastic_vector(x = model$mosquito$Y[, i] - BQ_success, prob = t(Psi_fail))
    model$mosquito$Y[, i] <- BQ_move
  }

  # for M we need to take into account Y0 (newly infecteds)
  M_success <- rbinom(n = l+p, size = model$mosquito$M, prob = c(psiB, psiQ))
  M_fail <- model$mosquito$M - M_success
  Y0 <- rbinom(n = p, size = M_success[1:p], prob = kappa)
  M_success[1:p] <- M_success[1:p] - Y0

  M_move <- sample_stochastic_vector(x = M_success, prob = t(Psi_success))
  M_move <- M_move + sample_stochastic_vector(x = M_fail, prob = t(Psi_fail))
  model$mosquito$M <- M_move

  Y0 <- sample_stochastic_vector(x = Y0, prob = t(model$mosquito$Psi_bq))

  # update incubating mosquitoes and add newly infecteds
  model$mosquito$Y <- model$mosquito$Y %*% model$mosquito$EIP_shift
  model$mosquito$Y[(p+1):(l+p), EIP+1L] <- Y0

  # newly emerging adults
  lambda <- compute_emergents(model)
  model$mosquito$M[1:p] <- model$mosquito$M[1:p] + model$mosquito$Psi_qb %*% lambda

}


# compute values for blood feeding

#' @title Compute mosquito feeding rate for BQ model (\eqn{f})
#' @description Blood feeding rates are modeled as a Holling type 2 (rational) function of blood host availability.
#' \deqn{f(B) = f_x \frac{s_f B}{1+s_f B}}
#' Here \eqn{f_x} is the maximum blood feeding rate and \eqn{s_f} is a scaling parameter.
#' @inheritParams compute_f
#' @return a vector of length `p` giving the per-capita blood feeding rate of mosquitoes in each blood feeding haunt
#' @export
compute_f.BQ <- function(model, B) {
  fx <- model$mosquito$fx
  sf <- model$mosquito$sf
  return(fx * ((sf * B) / (1 + (sf * B))))
}

#' @title Compute human blood feeding fraction for BQ model (\eqn{q})
#' @description The human blood feeding fraction is simply the proportion of human hosts.
#' @inheritParams compute_q
#' @return a vector of length `p` giving the proportion of bites taken on human hosts in each blood feeding haunt
#' @export
compute_q.BQ <- function(model, W, Wd, B) {
  return((W + Wd) / B)
}


#' @title Compute density of infective mosquitoes for BQ model (\eqn{Z})
#' @description This method returns `Z`.
#' @inheritParams compute_Z
#' @return a vector of length `p` giving the density of infected and infectious mosquitoes in each blood feeding haunt
#' @export
compute_Z.BQ <- function(model) {
  p <- model$global$p
  return(model$mosquito$Y[1:p, 1L])
}


# compute values for aquatic model

#' @title Compute number of eggs laid from oviposition for each aquatic habitat for BQ model
#' @description This method returns a vector of length `l`.
#' @inheritParams compute_oviposit
#' @return a vector of length `l` giving the total number of eggs laid by adult mosquitoes in each aquatic habitat
#' @details see [MicroMoB::compute_oviposit.BQ_deterministic] and [MicroMoB::compute_oviposit.BQ_stochastic]
#' @export
compute_oviposit.BQ <- function(model) {
  NextMethod()
}


#' @title Compute number of eggs laid from oviposition for each patch for deterministic RM model
#' @inheritParams compute_oviposit
#' @return a vector of length `l` giving the total number of eggs laid by adult mosquitoes in each aquatic habitat
#' @export
compute_oviposit.BQ_deterministic <- function(model) {
  l <- model$global$l
  p <- model$global$p
  psiQ <- model$mosquito$psiQ_mat[, tnow]
  Q <- model$mosquito$M[(p+1):(l+p), 1] + rowSums(model$mosquito$Y)[(p+1):(l+p)]
  return(model$mosquito$nu * psiQ * Q)
}


#' @title Compute number of eggs laid from oviposition for each patch for stochastic RM model
#' @inheritParams compute_oviposit
#' @return a vector of length `l` giving the total number of eggs laid by adult mosquitoes in each aquatic habitat
#' @importFrom stats rpois
#' @export
compute_oviposit.BQ_stochastic <- function(model) {
  l <- model$global$l
  p <- model$global$p
  psiQ <- model$mosquito$psiQ_mat[, tnow]
  Q <- model$mosquito$M[(p+1):(l+p), 1] + rowSums(model$mosquito$Y)[(p+1):(l+p)]
  return(rpois(n = p, lambda = model$mosquito$nu * psiQ * Q))
}
