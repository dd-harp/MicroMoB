#' @title Setup humans with SIR infection model
#' @description A simple SIR (Susceptible-Infected-Recovered) model
#' @param stochastic should the model update deterministically or stochastically?
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @param theta a time spent matrix
#' @param wf biting weights
#' @param H vector of strata population sizes
#' @param SIR a matrix giving S, I, R counts (columns) for each strata (rows)
#' @param b transmission efficiency (mosquito to human)
#' @param c transmission efficiency (human to mosquito)
#' @param gamma rate of recovery
#' @export
setup_humans_SIR <- function(model, stochastic, theta, wf = NULL, H, SIR, b = 0.55, c = 0.15, gamma = 1/5) {
  stopifnot(inherits(model, "MicroMoB"))
  stopifnot(inherits(theta, "matrix"))
  stopifnot(inherits(SIR, "matrix"))

  stopifnot(nrow(theta) >= ncol(theta))
  stopifnot(approx_equal(rowSums(theta), 1))

  n <- nrow(theta)
  p <- ncol(theta)
  stopifnot(p == model$global$p)

  stopifnot(length(H) == n)

  stopifnot(nrow(SIR) == n)
  stopifnot(ncol(SIR) == 3L)
  stopifnot(rowSums(SIR) == H)

  stopifnot(is.finite(c(b, c, gamma)))
  stopifnot(c(b, c, gamma) >= 0)
  stopifnot(c(b, c) <= 1)

  model$global$n <- n

  if (is.null(wf)) {
    wf <- rep(1, n)
  }

  stopifnot(length(wf) == n)
  stopifnot(is.finite(wf))
  stopifnot(wf >= 0)

  if (is.null(colnames(SIR))) {
    colnames(SIR) <- c("S", "I", "R")
  } else {
    stopifnot(colnames(SIR) == c("S", "I", "R"))
  }

  human_class <- c("SIR")
  if (stochastic) {
    human_class <- c(human_class, "SIR_stochastic")
  } else {
    human_class <- c(human_class, "SIR_deterministic")
  }

  model$human <- structure(list(), class = human_class)
  model$human$theta <- theta
  model$human$wf <- wf
  model$human$H <- H
  model$human$SIR <- SIR

  model$human$EIR <- rep(0, n)

  model$human$b <- b
  model$human$c <- c
  model$human$gamma <- gamma
}


# step (update)

#' @title Update SIR human model
#' @inheritParams step_humans
#' @export
step_humans.SIR <- function(model) {
  NextMethod()
}

#' @title Update SIR human model (deterministic)
#' @inheritParams step_humans
#' @importFrom stats pexp
#' @export
step_humans.SIR_deterministic <- function(model) {

  h <- model$human$EIR * model$human$b

  # compute differences: S
  S_leave <- model$human$SIR[, "S"] * pexp(q =  h)

  # compute differences: I
  I_leave <- model$human$SIR[, "I"] * pexp(q =  model$human$gamma)

  # update
  model$human$SIR[, "S"] <- model$human$SIR[, "S"] - S_leave
  model$human$SIR[, "I"] <- model$human$SIR[, "I"] + S_leave - I_leave
  model$human$SIR[, "R"] <- model$human$SIR[, "R"] + I_leave

  model$human$SIR <- pmax(model$human$SIR, 0)

}

#' @title Update SIR human model (stochastic)
#' @inheritParams step_humans
#' @importFrom stats pexp rbinom
#' @export
step_humans.SIR_stochastic <- function(model) {

  h <- model$human$EIR * model$human$b

  n <- model$global$n

  # compute differences: S
  S_leave <- rbinom(n = n, size = model$human$SIR[, "S"], prob = pexp(q =  h))

  # compute differences: I
  I_leave <- rbinom(n = n, size = model$human$SIR[, "I"], prob = pexp(q =  model$human$gamma))

  # update
  model$human$SIR[, "S"] <- model$human$SIR[, "S"] - S_leave
  model$human$SIR[, "I"] <- model$human$SIR[, "I"] + S_leave - I_leave
  model$human$SIR[, "R"] <- model$human$SIR[, "R"] + I_leave

}

#' @title Compute available humans for SIR model (\eqn{W})
#' @inheritParams compute_W
#' @export
compute_W.SIR <- function(model) {
  Psi <- model$human$theta
  W <- t(Psi) %*% (model$human$wf * model$human$H)
  return(as.vector(W))
}

#' @title Compute human biting weights for SIR model (\eqn{w_{f}})
#' @inheritParams compute_wf
#' @export
compute_wf.SIR <- function(model) {
  model$human$wf
}

#' @title Compute net infectiousness for SIR model (\eqn{x})
#' @inheritParams compute_x
#' @export
compute_x.SIR <- function(model) {
  X <- model$human$SIR[, "I"] / model$human$H
  return(as.vector(X * model$human$c))
}

#' @title Compute human population strata sizes for SIR model (\eqn{H})
#' @inheritParams compute_H
#' @export
compute_H.SIR <- function(model) {
  model$human$H
}


#' @title Compute time at risk matrix for SIR model (\eqn{\Psi})
#' @inheritParams compute_Psi
#' @export
compute_Psi.SIR <- function(model) {
  model$human$theta
}
