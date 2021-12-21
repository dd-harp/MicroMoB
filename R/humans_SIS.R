#' @title Setup humans with SIS pathogen model
#' @description A simple SIS model
#' @param stochastic should the model update deterministically or stochastically?
#' @param model an object from [MicroWNV::make_microWNV]
#' @param theta a time spent matrix
#' @param wf biting weights
#' @param H vector of strata population sizes
#' @param X number of infectious persons in each strata
#' @param b transmission efficiency (mosquito to human)
#' @param c transmission efficiency (human to mosquito)
#' @param r recovery rate (inverse of infectious duration)
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

  model$human$h <- rep(0, n)
  model$human$EIR <- rep(0, n)

  model$human$b <- b
  model$human$c <- c
  model$human$r <- r
}


# step (update)

#' @title Update SIS human model
#' @inheritParams step_humans
#' @export
step_humans.SIS <- function(model) {
  NextMethod()
}

#' @title Update SIS human model (deterministic)
#' @inheritParams step_humans
#' @importFrom stats pexp
#' @export
step_humans.SIS_deterministic <- function(model) {

  new_infections <- pexp(q = model$human$h) * (model$human$H - model$human$X)
  old_infections <- (1 - pexp(q = model$human$r)) * model$human$X

  model$human$X <- new_infections + old_infections

}

#' @title Update SIS human model (stochastic)
#' @inheritParams step_humans
#' @importFrom stats pexp rbinom
#' @export
step_humans.SIS_stochastic <- function(model) {

  n <- model$global$n
  new_infections <- rbinom(n = n, prob = pexp(q = model$human$h), size = model$human$H - model$human$X)
  old_infections <- rbinom(n = n, prob = 1 - pexp(q = model$human$r), size = model$human$X)

  model$human$X <- new_infections + old_infections

}


#' @title Compute available humans for SIS model
#' @inheritParams compute_W
#' @export
compute_W.SIS <- function(model) {
  Psi <- model$human$theta
  W <- t(Psi) %*% (model$human$wf * model$human$H)
  return(W)
}

#' @title Compute net infectiousness for SIR model
#' @inheritParams compute_x
#' @export
compute_x.SIS <- function(model) {
  X <- model$human$X / model$human$H
  return(X * model$human$c)
}

