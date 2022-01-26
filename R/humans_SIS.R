#' @title Setup humans with SIS pathogen model
#' @description A simple SIS (Susceptible-Infected-Susceptible) model
#' @param stochastic should the model update deterministically or stochastically?
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @param theta a time spent matrix
#' @param wf biting weights
#' @param H vector of strata population sizes
#' @param X number of infectious persons in each strata
#' @param b transmission efficiency (mosquito to human)
#' @param c transmission efficiency (human to mosquito)
#' @param r recovery rate (inverse of infectious duration)
#' @return no return value
#' @export
setup_humans_SIS <- function(model, stochastic, theta, wf = NULL, H, X, b = 0.55, c = 0.15, r = 1/200) {
  stopifnot(inherits(model, "MicroMoB"))
  stopifnot(inherits(theta, "matrix"))

  stopifnot(nrow(theta) >= ncol(theta))
  stopifnot(approx_equal(rowSums(theta), 1))

  n <- nrow(theta)
  p <- ncol(theta)
  stopifnot(p == model$global$p)

  stopifnot(length(H) == n)
  stopifnot(length(X) == n)
  stopifnot(X <= H)

  stopifnot(is.finite(c(b, c, r)))
  stopifnot(c(b, c, r) >= 0)
  stopifnot(c(b, c) <= 1)

  model$global$n <- n

  if (is.null(wf)) {
    wf <- rep(1, n)
  }

  stopifnot(length(wf) == n)
  stopifnot(is.finite(wf))
  stopifnot(wf >= 0)

  human_class <- c("SIS")
  if (stochastic) {
    human_class <- c(human_class, "SIS_stochastic")
  } else {
    human_class <- c(human_class, "SIS_deterministic")
  }

  model$human <- structure(list(), class = human_class)
  model$human$theta <- theta
  model$human$wf <- wf
  model$human$H <- H
  model$human$X <- X

  model$human$EIR <- rep(0, n)

  model$human$b <- b
  model$human$c <- c
  model$human$r <- r
}


# step (update)

#' @title Update SIS human model
#' @inheritParams step_humans
#' @return no return value
#' @export
step_humans.SIS <- function(model) {
  NextMethod()
}

#' @title Update SIS human model (deterministic)
#' @inheritParams step_humans
#' @return no return value
#' @importFrom stats pexp
#' @export
step_humans.SIS_deterministic <- function(model) {

  h <- model$human$EIR * model$human$b

  new_infections <- pexp(q = h) * (model$human$H - model$human$X)
  old_infections <- (1 - pexp(q = model$human$r)) * model$human$X

  model$human$X <- new_infections + old_infections

}

#' @title Update SIS human model (stochastic)
#' @inheritParams step_humans
#' @return no return value
#' @importFrom stats pexp rbinom
#' @export
step_humans.SIS_stochastic <- function(model) {

  h <- model$human$EIR * model$human$b
  n <- model$global$n

  new_infections <- rbinom(n = n, prob = pexp(q = h), size = model$human$H - model$human$X)
  old_infections <- rbinom(n = n, prob = 1 - pexp(q = model$human$r), size = model$human$X)

  model$human$X <- new_infections + old_infections

}


#' @title Compute available humans for SIS model (\eqn{W})
#' @inheritParams compute_W
#' @return a vector of length `p` giving the biting availability of human hosts at each patch
#' @export
compute_W.SIS <- function(model) {
  Psi <- model$human$theta
  W <- t(Psi) %*% (model$human$wf * model$human$H)
  return(as.vector(W))
}


#' @title Compute human biting weights for SIS model (\eqn{w_{f}})
#' @inheritParams compute_wf
#' @return a vector of length `n` giving the biting weights of human hosts in each stratum
#' @export
compute_wf.SIS <- function(model) {
  model$human$wf
}


#' @title Compute net infectiousness for SIS model (\eqn{x})
#' @inheritParams compute_x
#' @return a vector of length `n` giving the net infectiousness of human hosts in each stratum
#' @export
compute_x.SIS <- function(model) {
  x <- (model$human$X / model$human$H) * model$human$c
  return(as.vector(x))
}


#' @title Compute human population strata sizes for SIS model (\eqn{H})
#' @inheritParams compute_H
#' @return a vector of length `n` giving the size of each human population stratum
#' @export
compute_H.SIS <- function(model) {
  model$human$H
}


#' @title Compute time at risk matrix for SIS model (\eqn{\Psi})
#' @inheritParams compute_Psi
#' @return a matrix with `n` rows and `p` columns, the time at risk matrix
#' @export
compute_Psi.SIS <- function(model) {
  model$human$theta
}

