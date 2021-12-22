# classes and methods to implement a reasonably detailed RM model of mosquitoes

#' @title Setup generalized Ross-Macdonald mosquito model
#' @description This is a generalized RM model which allows for time varying EIP and
#' survival probability. It complies with the mosquito component interface, and may
#' be simulated deterministically or stochastically.
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @param stochastic should the model update deterministically or stochastically?
#' @param f the blood feeding rate
#' @param q the human blood feeding fraction
#' @param eip the Extrinsic Incubation Period, may either be a scalar, a vector of
#' length 365, or a vector of length equal to `tmax` in the `model` object from [MicroMoB::make_MicroMoB]
#' @param p daily survival probability, may either be a scalar, a vector of
#' length 365, or a vector of length equal to `tmax` in the `model` object from [MicroMoB::make_MicroMoB]
#' @param psi a mosquito dispersal matrix (rows must sum to 1)
#' @param M total mosquito density per patch (vector of length `p`)
#' @param Y density of incubating mosquitoes per patch (vector of length `p`)
#' @param Z density of infectious mosquitoes per patch (vector of length `p`)
#' @export
setup_mosquito_RM <- function(model, stochastic, f, q, eip, p, psi, M, Y, Z) {
  stopifnot(inherits(model, "MicroMoB"))
  stopifnot(is.logical(stochastic))

  if (stochastic) {
    M <- as.integer(M)
    Y <- as.integer(Y)
    Z <- as.integer(Z)
  }

  if (length(eip) == 1L) {
    stopifnot(is.finite(eip))
    stopifnot(eip > 0)
    eip_vec <- rep(eip, model$global$tmax)
  } else if(length(eip) == 365L) {
    stopifnot(is.finite(eip))
    stopifnot(eip > 0)
    ix <- (1:model$global$tmax) %% 365L
    ix[which(ix == 0L)] <- 365L
    eip_vec <- eip[ix]
  } else if(length(eip) == model$global$tmax) {
    stopifnot(is.finite(eip))
    stopifnot(eip > 0)
    eip_vec <- eip
  } else {
    stop("incorrect length of eip vector")
  }

  maxEIP <- max(eip_vec)

  if (length(p) == 1L) {
    stopifnot(is.finite(p))
    stopifnot(p > 0)
    p_vec <- rep(p, model$global$tmax)
  } else if(length(p) == 365L) {
    stopifnot(is.finite(p))
    stopifnot(p > 0)
    ix <- (1:model$global$tmax) %% 365L
    ix[which(ix == 0L)] <- 365L
    p_vec <- p[ix]
  } else if(length(p) == model$global$tmax) {
    stopifnot(is.finite(p))
    stopifnot(p > 0)
    p_vec <- p
  } else {
    stop("incorrect length of p vector")
  }

  stopifnot(p_vec <= 1)
  stopifnot(p_vec >= 0)
  stopifnot(!is.null(p_vec))

  stopifnot(dim(psi) == model$global$p)
  stopifnot(approx_equal(rowSums(psi), 1))

  mosy_class <- c("RM")
  if (stochastic) {
    mosy_class <- c(mosy_class, "RM_stochastic")
  } else {
    mosy_class <- c(mosy_class, "RM_deterministic")
  }

  if (length(f) == 1L) {
    f <- rep(f, model$global$p)
  } else {
    stopifnot(length(f) == model$global$p)
  }

  if (length(q) == 1L) {
    q <- rep(q, model$global$p)
  } else {
    stopifnot(length(q) == model$global$p)
  }

  model$mosquito <- structure(list(), class = mosy_class)
  model$mosquito$f <- f
  model$mosquito$q <- q
  model$mosquito$eip <- eip_vec
  model$mosquito$maxEIP <- maxEIP
  model$mosquito$p <- p_vec
  model$mosquito$psi <- psi

  model$mosquito$kappa <- rep(0, model$global$p)

  stopifnot(length(M) == model$global$p)
  stopifnot(length(Y) == model$global$p)
  stopifnot(length(Z) == model$global$p)
  stopifnot(is.finite(M))
  stopifnot(is.finite(Y))
  stopifnot(is.finite(Z))
  stopifnot(M >= Y)
  stopifnot(Y >= Z)

  model$mosquito$M <- M # mosquito density
  model$mosquito$Y <- Y # infected (incubating)
  model$mosquito$Z <- Z # infectious
  model$mosquito$ZZ <- matrix(data = 0, nrow = maxEIP, ncol = model$global$p) # each row is the number that will be added to the infectious state on that day
  if (any(Y > 0)) {
    if (stochastic) {
      ZZ <- Y - Z
      model$mosquito$ZZ[1:maxEIP, ] <- vapply(X = ZZ, FUN = function(z) {distribute(n = z, p = maxEIP)}, FUN.VALUE = numeric(maxEIP))
    } else {
      model$mosquito$ZZ[1:maxEIP, ] <- as.matrix(replicate(n = maxEIP, expr = (Y - Z)/maxEIP))
    }
  }

  ZZ_shift <- matrix(0, nrow = maxEIP, ncol = maxEIP)
  ZZ_shift[1:(maxEIP-1), 2:maxEIP] <- diag(maxEIP-1)
  model$mosquito$ZZ_shift <- ZZ_shift

}


# update mosquitoes over one time step

#' @title Update Ross-Macdonald mosquitoes
#' @description This function dispatches on the second argument of `model$mosquito`
#' for stochastic or deterministic behavior.
#' @inheritParams step_mosquitoes
#' @details see [MicroMoB::step_mosquitoes.RM_deterministic] and [MicroMoB::step_mosquitoes.RM_stochastic]
#' @export
step_mosquitoes.RM <- function(model) {
  NextMethod()
}

#' @title Update Ross-Macdonald mosquitoes (deterministic)
#' @inheritParams step_mosquitoes
#' @export
step_mosquitoes.RM_deterministic <- function(model) {

  # parameters
  tnow <- model$global$tnow
  EIP <- model$mosquito$eip[tnow]
  maxEIP <- model$mosquito$maxEIP
  p <- model$mosquito$p[tnow]
  psi <- model$mosquito$psi

  # newly infected mosquitoes
  a <- model$mosquito$f * model$mosquito$q
  Y0 <- a * model$mosquito$kappa * (model$mosquito$M - model$mosquito$Y)
  Y0 <- pmax(Y0, 0)

  # newly emerging adults
  lambda <- compute_emergents(model)

  # survival
  model$mosquito$M <- (p * model$mosquito$M) + lambda
  model$mosquito$Y <- p * (model$mosquito$Y + Y0)
  model$mosquito$Z <- p * model$mosquito$Z
  model$mosquito$ZZ <- p * model$mosquito$ZZ

  # dispersal
  model$mosquito$M <- model$mosquito$M %*% psi
  model$mosquito$Y <- model$mosquito$Y %*% psi
  model$mosquito$Z <- (model$mosquito$Z + model$mosquito$ZZ[1, ]) %*% psi
  model$mosquito$ZZ <- model$mosquito$ZZ %*% psi

  # ZZ[t, ] is the number of mosquitoes that become infectious in each patch t days from now.
  model$mosquito$ZZ <- model$mosquito$ZZ_shift %*% model$mosquito$ZZ
  model$mosquito$ZZ[EIP, ] <- model$mosquito$ZZ[EIP, ] + ((Y0 * p) %*% psi)

  # make vectors
  model$mosquito$M <- as.vector(model$mosquito$M)
  model$mosquito$Y <- as.vector(model$mosquito$Y)
  model$mosquito$Z <- as.vector(model$mosquito$Z)

}

#' @title Update Ross-Macdonald mosquitoes (stochastic)
#' @inheritParams step_mosquitoes
#' @importFrom stats rbinom rmultinom
#' @importFrom extraDistr rmvhyper
#' @export
step_mosquitoes.RM_stochastic <- function(model) {

  # parameters
  tnow <- model$global$tnow
  EIP <- model$mosquito$eip[tnow]
  maxEIP <- model$mosquito$maxEIP
  p <- model$mosquito$p[tnow]
  psi <- model$mosquito$psi
  n_patch <- model$global$p

  # newly emerging adults
  lambda <- compute_emergents(model)
  lambda <- sample_stochastic_vector(x = lambda, prob = psi)

  # newly infected mosquitoes
  a <- model$mosquito$f * model$mosquito$q
  model$mosquito$Y <- model$mosquito$Z + colSums(model$mosquito$ZZ)
  h <- a * model$mosquito$kappa
  Y0 <- rbinom(n = n_patch, size = model$mosquito$M - model$mosquito$Y, prob = h)

  # survival
  z_deaths <- rbinom(n = n_patch, size = model$mosquito$Z, prob = 1 - p)
  zz_deaths <- rbinom(n = n_patch, size = colSums(model$mosquito$ZZ), prob = 1 - p)
  y0_deaths <- rbinom(n = n_patch, size = Y0, prob = 1 - p)
  s_deaths <- rbinom(n = n_patch, size = model$mosquito$M - model$mosquito$Y - Y0, prob = 1 - p) # S (susceptible) = M - Y - Y0

  # adjust zz death to be bin (day) specific
  zz_deaths <- vapply(X = 1:n_patch, FUN = function(x) {
    if (zz_deaths[x] > 0) {
      # scatter deaths across incubating bins (equiprobable for each bin)
      rmvhyper(nn = 1, n = model$mosquito$ZZ[, x], k = zz_deaths[x])
    } else {
      rep(0, maxEIP)
    }
  }, FUN.VALUE = numeric(maxEIP), USE.NAMES = FALSE)

  model$mosquito$ZZ <- model$mosquito$ZZ - zz_deaths
  model$mosquito$Z <- model$mosquito$Z - z_deaths
  Y0 <- Y0 - y0_deaths
  S <- (model$mosquito$M - model$mosquito$Y - Y0) - s_deaths

  # dispersal
  model$mosquito$Z <- sample_stochastic_vector(x = model$mosquito$Z, prob = psi)
  model$mosquito$ZZ <- sample_stochastic_matrix(x = model$mosquito$ZZ, prob = psi)
  Y0 <- sample_stochastic_vector(x = Y0, prob = psi)
  S <- sample_stochastic_vector(x = S, prob = psi)

  # add newly incubating
  model$mosquito$Y <- colSums(model$mosquito$ZZ) + model$mosquito$Z + Y0

  # add newly infectious
  model$mosquito$Z <- model$mosquito$Z + model$mosquito$ZZ[1, ]

  # add newly emerging
  model$mosquito$M <- model$mosquito$Y + S + lambda

  # ZZ[t, ] is the number of mosquitoes that become infectious in each patch t days from now.
  model$mosquito$ZZ <- model$mosquito$ZZ_shift %*% model$mosquito$ZZ
  model$mosquito$ZZ[EIP, ] <- model$mosquito$ZZ[EIP, ] + Y0

}


# compute values for blood feeding

#' @title Compute mosquito feeding rate for RM model (\eqn{f})
#' @description This method simply returns the `f` parameter of the mosquito object,
#' because the RM model assumes a constant blood feeding rate.
#' @inheritParams compute_f
#' @export
compute_f.RM <- function(model, B) {
  model$mosquito$f
}

#' @title Compute human blood feeding fraction for RM model (\eqn{q})
#' @description This method simply returns the `q` parameter of the mosquito object,
#' because the RM model assumes a constant fraction of blood meals are taken on
#' human hosts.
#' @inheritParams compute_q
#' @export
compute_q.RM <- function(model, W, Wd, B) {
  model$mosquito$q
}


#' @title Compute density of infective mosquitoes for RM model (\eqn{Z})
#' @description This method returns `Z`.
#' @inheritParams compute_Z
#' @export
compute_Z.RM <- function(model) {
  model$mosquito$Z
}

