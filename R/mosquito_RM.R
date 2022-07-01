# classes and methods to implement a reasonably detailed RM model of mosquitoes

#' @title Setup generalized Ross-Macdonald mosquito model
#' @description This is a generalized RM model which allows for time varying EIP and
#' survival probability. It complies with the mosquito component interface, and may
#' be simulated deterministically or stochastically.
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @param stochastic should the model update deterministically or stochastically?
#' @param f the blood feeding rate
#' @param q the human blood feeding fraction
#' @param eip the Extrinsic Incubation Period (may be time varying see [MicroMoB::time_varying_parameter])
#' @param p daily survival probability (may be time and patch varying see [MicroMoB::time_patch_varying_parameter])
#' @param psi a mosquito dispersal matrix (rows must sum to 1)
#' @param nu number of eggs laid per oviposition
#' @param M total mosquito density per patch (vector of length `p`)
#' @param Y density of incubating mosquitoes per patch (vector of length `p`)
#' @param Z density of infectious mosquitoes per patch (vector of length `p`)
#' @param N `l` by `p` matrix describing how eggs from mosquitoes in patches are
#' distributed amongst aquatic habitats. If `NULL` it is the identity matrix of dimension
#' `l`.
#' @return no return value
#' @importFrom stats rmultinom
#' @export
setup_mosquito_RM <- function(model, stochastic, f = 0.3, q = 0.9, eip, p, psi, nu = 25, M, Y, Z, N = NULL) {
  stopifnot(inherits(model, "MicroMoB"))
  stopifnot(is.logical(stochastic))

  tmax <- model$global$tmax

  if (is.null(N)) {
    # because in general we don't know how to map eggs in patches to habitats,
    # we do not allow l != p if N is not provided
    stopifnot(model$global$l == model$global$p)
    N <- diag(model$global$l)
  } else {
    stopifnot(nrow(N) == model$global$l)
    stopifnot(ncol(N) == model$global$p)
    stopifnot(approx_equal(colSums(N), 1))
  }

  if (stochastic) {
    M <- as.integer(M)
    Y <- as.integer(Y)
    Z <- as.integer(Z)
  }

  eip_vec <- time_varying_parameter(param = eip, tmax = tmax)

  maxEIP <- max(eip_vec)

  p_vec <- time_patch_varying_parameter(param = p, p = model$global$p, tmax = tmax)

  stopifnot(p_vec <= 1)
  stopifnot(p_vec >= 0)
  stopifnot(!is.null(p_vec))

  stopifnot(dim(psi) == model$global$p)
  stopifnot(approx_equal(rowSums(psi), 1))

  stopifnot(length(f) == 1L)
  stopifnot(length(q) == 1L)
  stopifnot(length(nu) == 1L)

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
  model$mosquito$nu <- nu
  model$mosquito$N <- N # oviposition matrix

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

    ZZ <- Y - Z
    p_surv <- model$mosquito$p[, 1L]
    denom <- vapply(X = 0:(maxEIP-1), FUN = function(n) {p_surv^n}, FUN.VALUE = numeric(model$global$p))
    num <- vapply(X = (maxEIP-1):0, FUN = function(n) {ZZ*(p_surv^n)}, FUN.VALUE = numeric(model$global$p))
    if (model$global$p > 1) {
      denom <- rowSums(denom)
      num <- t(num)
    } else {
      denom <- sum(denom)
      num <- as.matrix(num)
    }

    if (stochastic) {
      for (i in 1:model$global$p) {
        if (ZZ[i] > 0) {
          probs <- num[, i] / denom[i]
          model$mosquito$ZZ[, i] <- rmultinom(n = 1L, size = ZZ[i], prob = probs)
        }
      }
    } else {
      for (i in 1:model$global$p) {
        model$mosquito$ZZ[, i] <- num[, i] / denom[i]
      }
    }
  }

  ZZ_shift <- matrix(0, nrow = maxEIP, ncol = maxEIP)
  ZZ_shift[1:(maxEIP-1), 2:maxEIP] <- diag(maxEIP-1)
  model$mosquito$ZZ_shift <- ZZ_shift

}


#' @title Get parameters for generalized Ross-Macdonald mosquito model
#' @description The JSON config file should have 8 entries:
#'  * stochastic: a boolean value
#'  * f: scalar
#'  * q: scalar
#'  * eip: scalar or vector; see [MicroMoB::time_varying_parameter] for valid formats
#'  * p: scalar or vector; see [MicroMoB::time_varying_parameter] for valid formats
#'  * psi: matrix
#'  * nu: scalar
#'  * M: vector
#'  * Y: vector
#'  * Z: vector
#'
#' For interpretation of the entries, please read [MicroMoB::setup_mosquito_RM].
#' @param path a file path to a JSON file
#' @return a named [list]
#' @importFrom jsonlite read_json
#' @examples
#' # to see an example of proper JSON input, run the following
#' library(jsonlite)
#' t <- 10 # days to simulate
#' p <- 5 # number of patches
#' EIP <-  rep(5, t)
#' p_surv <- 0.95
#' psi <- matrix(rexp(p^2), nrow = p, ncol = p)
#' psi <- psi / rowSums(psi)
#' par <- list(
#'  "stochastic" = FALSE,
#'  "f" = 0.3,
#'  "q" = 0.9,
#'  "eip" = EIP,
#'  "p" = p_surv,
#'  "psi" = psi,
#'  "nu" = 20,
#'  "M" = rep(100, p),
#'  "Y" = rep(20, p),
#'  "Z" = rep(5, p)
#' )
#' toJSON(par, pretty = TRUE)
#' @export
get_config_mosquito_RM <- function(path) {
  pars <- read_json(path = file.path(path), simplifyVector = TRUE)

  stopifnot(length(pars) == 10L)
  stopifnot(is.logical(pars$stochastic))

  stopifnot(is.numeric(pars$f))
  stopifnot(is.numeric(pars$q))

  stopifnot(is.numeric(pars$eip))
  stopifnot(is.vector(pars$eip))

  stopifnot(is.numeric(pars$p))
  stopifnot(is.vector(pars$p))

  stopifnot(is.numeric(pars$psi))
  stopifnot(is.matrix(pars$psi))

  stopifnot(is.numeric(pars$nu))

  stopifnot(is.numeric(pars$M))
  stopifnot(is.numeric(pars$Y))
  stopifnot(is.numeric(pars$Z))
  stopifnot(is.vector(pars$M))
  stopifnot(is.vector(pars$Y))
  stopifnot(is.vector(pars$Z))

  return(pars)
}


# output

#' @title Get output for Ross-Macdonald mosquito populations
#' @description Return a [data.frame].
#' @inheritParams output_mosquitoes
#' @return a [data.frame] with columns `M` (all adult mosquitoes), `Y` (infected mosquitoes), and `Z` (infectious mosquitoes), and rows
#' correspond to places.
#' @export
output_mosquitoes.RM <- function(model) {
  data.frame(M = model$mosquito$M, Y = model$mosquito$Y, Z = model$mosquito$Z)
}


# update mosquitoes over one time step

#' @title Update Ross-Macdonald mosquitoes
#' @description This function dispatches on the second argument of `model$mosquito`
#' for stochastic or deterministic behavior.
#' @inheritParams step_mosquitoes
#' @return no return value
#' @details see [MicroMoB::step_mosquitoes.RM_deterministic] and [MicroMoB::step_mosquitoes.RM_stochastic]
#' @export
step_mosquitoes.RM <- function(model) {
  NextMethod()
}

#' @title Update Ross-Macdonald mosquitoes (deterministic)
#' @inheritParams step_mosquitoes
#' @return no return value
#' @export
step_mosquitoes.RM_deterministic <- function(model) {

  # parameters
  tnow <- model$global$tnow
  EIP <- model$mosquito$eip[tnow]
  p <- model$mosquito$p[, tnow]
  psi <- model$mosquito$psi

  # newly infected mosquitoes
  a <- model$mosquito$f * model$mosquito$q
  Y0 <- a * model$mosquito$kappa * (model$mosquito$M - model$mosquito$Y)
  Y0 <- pmax(Y0, 0)

  # survival
  model$mosquito$M <- p * model$mosquito$M
  model$mosquito$Y <- p * (model$mosquito$Y + Y0)
  model$mosquito$Z <- p * model$mosquito$Z
  model$mosquito$ZZ <- (matrix(data = 1, nrow = nrow(model$mosquito$ZZ), ncol = 1) %*% p) * model$mosquito$ZZ

  # dispersal
  model$mosquito$M <- model$mosquito$M %*% psi
  model$mosquito$Y <- model$mosquito$Y %*% psi
  model$mosquito$Z <- (model$mosquito$Z + model$mosquito$ZZ[1, ]) %*% psi
  model$mosquito$ZZ <- model$mosquito$ZZ %*% psi

  # ZZ[t, ] is the number of mosquitoes that become infectious in each patch t days from now.
  model$mosquito$ZZ <- model$mosquito$ZZ_shift %*% model$mosquito$ZZ
  model$mosquito$ZZ[EIP, ] <- model$mosquito$ZZ[EIP, ] + ((Y0 * p) %*% psi)

  # newly emerging adults
  lambda <- compute_emergents(model)

  # make vectors
  model$mosquito$M <- as.vector(model$mosquito$M) + lambda
  model$mosquito$Y <- as.vector(model$mosquito$Y)
  model$mosquito$Z <- as.vector(model$mosquito$Z)

}


#' @title Update Ross-Macdonald mosquitoes (stochastic)
#' @inheritParams step_mosquitoes
#' @return no return value
#' @importFrom stats rbinom
#' @export
step_mosquitoes.RM_stochastic <- function(model) {

  # parameters
  tnow <- model$global$tnow
  EIP <- model$mosquito$eip[tnow]
  p <- model$mosquito$p[, tnow]
  psi <- model$mosquito$psi
  n_patch <- model$global$p
  maxEIP <- model$mosquito$maxEIP

  # newly infected mosquitoes
  a <- model$mosquito$f * model$mosquito$q
  Y0 <- rbinom(n = n_patch, size = model$mosquito$M - model$mosquito$Y, prob = a * model$mosquito$kappa)

  # susceptible mosquitoes
  S <- model$mosquito$M - model$mosquito$Z - colSums(model$mosquito$ZZ) - Y0

  # survival
  S <- rbinom(n = n_patch, size = S, prob = p)
  Y0 <- rbinom(n = n_patch, size = Y0, prob = p)
  model$mosquito$ZZ <- matrix(
    data = rbinom(n = prod(dim(model$mosquito$ZZ)), size = as.vector(model$mosquito$ZZ), prob = p),
    nrow = maxEIP, ncol = n_patch
  )
  model$mosquito$Z <- rbinom(n = n_patch, size = model$mosquito$Z, prob = p)

  # dispersal
  S <- sample_stochastic_vector(x = S, prob = psi)
  Y0 <- sample_stochastic_vector(x = Y0, prob = psi)
  model$mosquito$ZZ <- sample_stochastic_matrix(x = model$mosquito$ZZ, prob = psi)
  model$mosquito$Z <- sample_stochastic_vector(x = model$mosquito$Z, prob = psi)

  # add newly incubating
  model$mosquito$Z <- model$mosquito$Z + model$mosquito$ZZ[1L, ]

  # ZZ[t, ] is the number of mosquitoes that become infectious in each patch t days from now.
  model$mosquito$ZZ <- model$mosquito$ZZ_shift %*% model$mosquito$ZZ
  model$mosquito$ZZ[EIP, ] <- model$mosquito$ZZ[EIP, ] + Y0

  # newly emerging adults
  lambda <- compute_emergents(model)

  # make vectors
  model$mosquito$Y <- colSums(model$mosquito$ZZ) + model$mosquito$Z
  model$mosquito$M <- model$mosquito$Y + S + lambda

}


# compute values for blood feeding

#' @title Compute mosquito feeding rate for RM model (\eqn{f})
#' @description This method simply returns the `f` parameter of the mosquito object,
#' because the RM model assumes a constant blood feeding rate.
#' @inheritParams compute_f
#' @return a vector of length `p` giving the per-capita blood feeding rate of mosquitoes in each patch
#' @export
compute_f.RM <- function(model, B) {
  model$mosquito$f
}

#' @title Compute human blood feeding fraction for RM model (\eqn{q})
#' @description This method simply returns the `q` parameter of the mosquito object,
#' because the RM model assumes a constant fraction of blood meals are taken on
#' human hosts.
#' @inheritParams compute_q
#' @return a vector of length `p` giving the proportion of bites taken on human hosts in each patch
#' @export
compute_q.RM <- function(model, W, Wd, B) {
  model$mosquito$q
}


#' @title Compute density of infective mosquitoes for RM model (\eqn{Z})
#' @description This method returns `Z`.
#' @inheritParams compute_Z
#' @return a vector of length `p` giving the density of infected and infectious mosquitoes in each patch
#' @export
compute_Z.RM <- function(model) {
  model$mosquito$Z
}


# compute values for aquatic model

#' @title Compute number of eggs laid from oviposition for each patch for RM model
#' @description This method returns a vector of length `p`.
#' @inheritParams compute_oviposit
#' @return a vector of length `p` giving the total number of eggs laid by adult mosquitoes in each patch
#' @details see [MicroMoB::compute_oviposit.RM_deterministic] and [MicroMoB::compute_oviposit.RM_stochastic]
#' @export
compute_oviposit.RM <- function(model) {
  NextMethod()
}


#' @title Compute number of eggs laid from oviposition for each patch for deterministic RM model
#' @inheritParams compute_oviposit
#' @return a vector of length `p` giving the total number of eggs laid by adult mosquitoes in each patch
#' @export
compute_oviposit.RM_deterministic <- function(model) {
  return(as.vector(model$mosquito$N %*% (model$mosquito$nu * model$mosquito$f * model$mosquito$M)))
}


#' @title Compute number of eggs laid from oviposition for each patch for stochastic RM model
#' @inheritParams compute_oviposit
#' @return a vector of length `l` giving the total number of eggs laid by adult mosquitoes in each patch
#' @importFrom stats rpois
#' @export
compute_oviposit.RM_stochastic <- function(model) {
  eggs <- rpois(n = model$global$p, lambda = model$mosquito$nu * model$mosquito$f * model$mosquito$M)
  return(as.vector(model$mosquito$N %*% eggs))
}
