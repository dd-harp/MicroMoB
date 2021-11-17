# time spent models (dt: day divided into d equal fractions, and day: just one day)

#' @title Setup a time spent model
#' @description Setup a time spent model. The model object must have already
#' been initialized with a human object (see [MicroMoB::setup.human]). This adds
#' a [list] to `model` named "tisp" (Time spent).
#' @section Daily time spent model: see [MicroMoB::setup.timespent.day]
#' @section Fractional day time spent model: see [MicroMoB::setup.timespent.dt]
#' @param type a character in `c("day")`
#' @param model a model object (an [environment])
#' @param ... other arguments to be passed to type methods
#' @export
setup.timespent <- function(type, model, ...) {
  stopifnot(inherits(model, "environment"))
  stopifnot(!is.null(model$human))
  stopifnot(type %in% c("day", "dt"))
  tisp <- structure(list(), class = type)
  UseMethod("setup.timespent", tisp)
}

#' @title Setup daily time spent model
#' @description A model of where people spend their time over the day.
#' @inheritParams setup.timespent
#' @param theta a time spent matrix, if `NULL` the identity matrix is used
#' (everyone stays at their home patch), otherwise `theta` should be a
#' matrix whose rows sum to \eqn{\le 1} giving the daily time spent distribution.
#' @export
setup.timespent.day <- function(type, model, theta = NULL, ...) {

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
  tisp$Theta_t <- t(theta)
  model$tisp <- tisp
}

#' @title Setup fractional day time spent model
#' @inheritParams setup.timespent
#' @param theta a time spent matrix, if `NULL` the identity matrix is used
#' (everyone stays at their home patch)
#' @export
setup.timespent.dt <- function(type, model, theta, ...) {

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
  tisp$d <- d
  tisp$Theta_t <- array(data = do.call(c, lapply(X = theta, FUN = function(x){t(x)})), dim = c(p, n, d))
  model$tisp <- tisp
}


# compute time at risk (Psi)

#' @title Compute time at risk (\eqn{\Psi})
#' @param tisp a time spent object (see [MicroMoB::setup.timespent])
#' @param xi probabilities to initiate a feeding search during each fraction
#' of the day (must sum to 1)
#' @param t time (day)
#' @seealso [MicroMoB::compute_Psi.timespent.day] [MicroMoB::compute_Psi.timespent.dt]
#' @return a (transposed) time at risk [matrix] or [array], which should be attached to the `tisp`
#' object and stored for later computation (human availability and the biting distribution matrix)
#' @export
compute_Psi.timespent <- function(tisp, xi, t) {
  stopifnot(is.finite(t))
  stopifnot(sum(xi) == 1)
  UseMethod("compute_Psi.timespent", tisp)
}

#' @title Compute time at risk for daily time spent (\eqn{\Psi^{T}})
#' @inheritParams compute_Psi.timespent
#' @export
compute_Psi.timespent.day <- function(tisp, xi, t) {
  stopifnot(length(xi) == 1L)
  return(tisp$Theta_t * xi[1])
}

#' @title Compute time at risk for fractional daily time spent (\eqn{\Psi^{T}})
#' @inheritParams compute_Psi.timespent
#' @export
compute_Psi.timespent.dt <- function(tisp, xi, t) {
  stopifnot(length(xi) == tisp$d)
  Psi_t <- tisp$Theta_t * 0
  for (k in 1:tisp$d) {
    Psi_t[, , k] <- tisp$Theta_t[, , k] * xi[k]
  }
  return(Psi_t)
}


# compute human availability (W)

#' @title Compute human availability (W)
#' @param tisp an object from [MicroMoB::setup.timespent]
#' @param biteweight an object from [MicroMoB::setup.biteweight]
#' @param human an object from [MicroMoB::setup.human]
#' @param t time
#' @seealso [MicroMoB::compute_W.timespent.day] [MicroMoB::compute_W.timespent.dt]
#' @export
compute_W.timespent <- function(tisp, biteweight, human, t) {
  stopifnot(is.finite(t))
  UseMethod("compute_W.timespent", tisp)
}

#' @title Compute human availability for fractional daily time spent (W)
#' @description For each time period the human availability is \eqn{W\[t\] = \Psi\[t\]^{T} \cdot w_f H},
#' and the sum of them over the day is the overall human availability.
#' @inheritParams compute_W.timespent
#' @export
compute_W.timespent.dt <- function(tisp, biteweight, human, t) {
  wf <- compute.biteweight(biteweight = biteweight, t = t)
  W <- lapply(X = 1:tisp$d, FUN = function(k){
    tisp$Psi_t[, , k] %*% (wf * human$H)
  })
  W <- do.call("+", W)
  return(W)
}

#' @title Compute human availability for daily time spent (W)
#' @description The human availability is \eqn{\eqn{W = \Psi^{T} \cdot w_f H},}
#' @inheritParams compute_W.timespent
#' @export
compute_W.timespent.day<- function(tisp, biteweight, human, t) {
  wf <- compute.biteweight(biteweight = biteweight, t = t)
  W <- tisp$Psi_t %*% (wf * human$H)
  return(W)
}
