# time spent models (dt: day divided into d equal fractions, and day: just one day)

#' @title Setup a time spent model
#' @description Setup a time spent model. The model object must have already
#' been initialized with a human object (see [MicroMoB::setup_human]). This adds
#' a [list] to `model$human` named "timespent" (Time spent).
#' @seealso [MicroMoB::setup_timespent.day] [MicroMoB::setup_timespent.dt]
#' @param type a character in `c("day")`
#' @param model a model object (an [environment])
#' @param ... other arguments to be passed to type methods
#' @export
setup_timespent <- function(type, model, ...) {
  stopifnot(inherits(model, "environment"))
  stopifnot(!is.null(model$human))
  timespent <- structure(list(), class = type)
  UseMethod("setup_timespent", timespent)
}

#' @title Setup daily time spent model
#' @description A model of where people spend their time over the day.
#' @inheritParams setup_timespent
#' @param theta a time spent matrix, if `NULL` the identity matrix is used
#' (everyone stays at their home patch), otherwise `theta` should be a
#' matrix whose rows sum to \eqn{\le 1} giving the daily time spent distribution.
#' @export
setup_timespent.day <- function(type, model, theta = NULL, ...) {

  if (is.null(theta)) {
    if (model$human$n == model$human$p) {
      theta <- diag(model$human$n)
    } else {
      theta <- t(model$human$J)
    }
  } else {
    stopifnot(is.finite(theta))
    stopifnot(rowSums(theta) <= 1)
    stopifnot(rowSums(theta) > 0)
  }

  p <- nrow(model$human$J)
  n <- length(model$human$H)
  stopifnot(nrow(theta) == n)
  stopifnot(ncol(theta) == p)

  # store transposed theta
  timespent$Theta_t <- t(theta)
  model$human$timespent <- timespent
}

#' @title Setup fractional day time spent model
#' @inheritParams setup_timespent
#' @param theta a time spent matrix, if `NULL` the identity matrix is used
#' (everyone stays at their home patch)
#' @export
setup_timespent.dt <- function(type, model, theta, ...) {

  stopifnot(inherits(theta, "list"))
  stopifnot(length(theta) > 1)

  p <- nrow(model$human$J)
  n <- length(model$human$H)
  d <- length(theta)

  for (k in 1:d) {
    stopifnot(nrow(theta[[k]]) == n)
    stopifnot(ncol(theta[[k]]) == p)
    stopifnot(is.finite(theta[[k]]))
    stopifnot(rowSums(theta[[k]]) <= 1)
    stopifnot(rowSums(theta[[k]]) > 0)
  }

  # store transposed theta
  timespent$d <- d
  timespent$Theta_t <- array(data = do.call(c, lapply(X = theta, FUN = function(x){t(x)})), dim = c(p, n, d))
  model$human$timespent <- timespent
}

#' @export
setup_timespent.default <- function(type, model, theta, ...) {
  stop("setup_timespent has no method for dispatch type ", type)
}


# compute time at risk (Psi)

#' @title Compute time at risk (\eqn{\Psi^{T}})
#' @description Dispatch is done on the class of `human$timespent`
#' @param human an object from [MicroMoB::setup_human]
#' @param xi probabilities to initiate a feeding search during each fraction
#' of the day (must sum to 1)
#' @param t time (day)
#' @seealso [MicroMoB::compute_Psi.day] [MicroMoB::compute_Psi.dt]
#' @return a (transposed) time at risk [matrix] or [array], which should be attached to the `timespent`
#' object and stored for later computation (human availability and the biting distribution matrix)
#' @export
compute_Psi <- function(human, xi, t) {
  stopifnot(is.finite(t))
  stopifnot(sum(xi) == 1)
  UseMethod("compute_Psi", human$timespent)
}

#' @title Compute time at risk for daily time spent (\eqn{\Psi^{T}})
#' @inheritParams compute_Psi
#' @return a matrix of dimension \eqn{p \times n}
#' @export
compute_Psi.day <- function(human, xi, t) {
  stopifnot(length(xi) == 1L)
  return(human$timespent$Theta_t * xi[1])
}

#' @title Compute time at risk for fractional daily time spent (\eqn{\Psi^{T}})
#' @inheritParams compute_Psi
#' @return an array of dimension \eqn{p \times n \time d}
#' @export
compute_Psi.dt <- function(human, xi, t) {
  stopifnot(length(xi) == human$timespent$d)
  Psi_t <- human$timespent$Theta_t * 0
  for (k in 1:human$timespent$d) {
    Psi_t[, , k] <- human$timespent$Theta_t[, , k] * xi[k]
  }
  return(Psi_t)
}

#' @export
compute_Psi.default <- function(human, xi, t) {
  stop("compute_Psi has no method for dispatch type ", class(human$timespent))
}


# compute human availability (W)

#' @title Compute human availability (W)
#' @param human an object from [MicroMoB::setup_human]
#' @param Psi_t either a matrix from [MicroMoB::compute_Psi.day] or an array
#' from [MicroMoB::compute_Psi.dt]
#' @param t time
#' @seealso [MicroMoB::compute_W.day] [MicroMoB::compute_W.dt]
#' @export
compute_W <- function(human, Psi_t, t) {
  stopifnot(is.finite(t))
  UseMethod("compute_W", human$timespent)
}

#' @title Compute human availability for fractional daily time spent (W)
#' @description For each time period the human availability is \eqn{W\[t\] = \Psi\[t\]^{T} \cdot w_f H},
#' and the sum of them over the day is the overall human availability.
#' @inheritParams compute_W
#' @return a matrix of dimension \eqn{p \times d}
#' @export
compute_W.dt <- function(human, Psi_t, t) {
  wf <- compute_wf(biteweight = human$biteweight, t = t)
  W <- vapply(X = 1:human$timespent$d, FUN = function(k){
    Psi_t[, , k] %*% (wf * human$H)
  }, FUN.VALUE = numeric(human$p))
  return(W)
}

#' @title Compute human availability for daily time spent (W)
#' @description The human availability is \eqn{\eqn{W = \Psi^{T} \cdot w_f H},}
#' @inheritParams compute_W
#' @return a vector of length \eqn{p}
#' @export
compute_W.day<- function(human, Psi_t, t) {
  wf <- compute_wf(biteweight = human$biteweight, t = t)
  W <- Psi_t %*% (wf * human$H)
  return(W)
}

#' @export
compute_W.default <- function(human, Psi_t, t) {
  stop("compute_W has no method for dispatch type ", class(human$timespent))
}
