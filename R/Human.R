

#' @title Helper function for lumped population strata
#' @description If input is given as a vector of population sizes per-strata, lumped
#' over patches, and a separate matrix whose columns describe how each strata is
#' distributed over patches, this function calculates the residency matrix and
#' population size for the overall stratification of both residency and strata.
#' @param H_strata a vector of population size by strata
#' @param J_strata a matrix chose columns sum to one giving the distribution of
#' strata populations over patches
#' @return a [list] with three elements:
#'  * `assignment_indices`: provides a mapping from patch (rows) and strata (columns)
#'     into the "unrolled" vector `H`
#'  * `J`: the residency matrix mapping elements in `H` to patches
#'  * `H`: the overall population distribution over strata and patches
#' @examples
#' # taken from package tests
#' J <- matrix(
#'    c(0.3, 0.5, 0.2,
#'    0.1, 0.6, 0.3), nrow = 3, ncol = 2, byrow = FALSE
#' )
#' H <- c(50, 60)
#' H_overall <- J %*% diag(H)
#' residency <- strata_to_residency(H_strata = H, J_strata = J)
#' @export
strata_to_residency <- function(H_strata, J_strata) {

  stopifnot(length(H_strata) == ncol(J_strata))
  stopifnot(colSums(J_strata) == 1)

  p <- nrow(J_strata)
  s <- ncol(J_strata)
  n <- p*s

  # we need to determine the actual residency matrix
  pop <- list()

  pop$assignment_indices <- matrix(data = 1:n, nrow = p, ncol = s, byrow = FALSE)

  pop$J <- matrix(0, nrow = p, ncol = n)
  for (j in 1:p) {
    pop$J[j, pop$assignment_indices[j, ]] <- 1
  }
  pop$H <- as.vector(J_strata %*% diag(H_strata))

  return(pop)
}

#' @title Setup a human object
#' @description Setup a human object. The `model` object will have a list named `human`
#' added to it.
#' @section strata: see [MicroMoB::setup.human.strata]
#' @param type a character in `c("simple", "strata")`
#' @param model a model object (an [environment])
#' @param ... other arguments to be passed to type methods
#' @export
setup.human <- function(type, model, ...) {
  stopifnot(inherits(model, "environment"))
  stopifnot(type %in% c("strata"))
  pop <- structure(list(), class = type)
  UseMethod("setup.human", pop)
}

#' @title Setup a human model with strata
#' @description This sets up a human model object.
#' @inheritParams setup.human
#' @param H a vector of human population sizes
#' @param J a matrix whose columns assign human strata to patches (rows); the
#' columns must all sum to one. If `J` is `NULL` then a diagonal matrix of ones
#' is used; this assumes that the only level of population stratification is
#' patch of residence.
#' @export
setup.human.strata <- function(type, model, H, J = NULL, ...) {

  stopifnot(length(H) > 0)
  stopifnot(is.finite(H))
  stopifnot(H >= 0)

  # if J is not provided, we have to assume that H one to one maps to patches
  if (is.null(J)) {

    p <- length(H)
    n <- p

    J <- diag(n)

  } else {

    p <- nrow(J)
    n <- length(H)

    stopifnot(is_binary(J))
    stopifnot(ncol(J) == n)
    stopifnot(colSums(J) == 1)
    stopifnot(n >= p)

  }

  pop$J <- J
  pop$H <- H

  pop$n <- n
  pop$p <- p

  model$human <- pop
}




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
  tisp$theta_t <- t(theta)
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
  tisp$theta_t <- array(data = do.call(c, lapply(X = theta, FUN = function(x){t(x)})), dim = c(p, n, d))
  model$tisp <- tisp
}


# computes W (human availability)



#' #' @title Compute human availability (W)
#' @description do thing
#' @section dispatching on a type of `timespent`
#' @param tisp an object from [MicroMoB::setup.timespent]
#' @param biteweight an object from [MicroMoB::setup.biteweight]
#' @param human an object from [MicroMoB::setup.human]
#' @param xi a vector of probabilities of mosquito feeding initiation
#' @param t time
#' @export
compute.timespent <- function(tisp, biteweight, human, xi, t) {
  stopifnot(is.finite(t))
  stopifnot(sum(xi) == 1)
  UseMethod("compute.timespent", tisp)
}

#' @rdname compute.timespent
#' @method compute.timespent dt
#' @export
compute.timespent.dt <- function(tisp, biteweight, human, xi, t) {
  stopifnot(length(xi) == tisp$d)
  wt <- compute.biteweight(biteweight = biteweight, t = t)
  W <- lapply(X = 1:tisp$d, FUN = function(k){
    # theta_t * xi is the TaR matrix
    tisp$theta_t[, , k] %*% (wt * human$H) * xi[k]
  })
  return(do.call("+", W))
}




#' @title Setup a time spent model
#' @description Setup a time spent model. The model object must have already
#' been initialized with a human object (see [MicroMoB::setup.human]). This adds
#' a [list] to `model` named "tisp" (Time spent).
#' @param type a character in `c("null")`
#' @param model a model object (an [environment])
#' @param ... other arguments to be passed to type methods
#' @export
setup.biteweight <- function(type, model, ...) {
  stopifnot(inherits(model, "environment"))
  stopifnot(!is.null(model$human))
  stopifnot(type %in% c("null"))
  biteweight <- structure(list(), class = type)
  UseMethod("setup.biteweight", biteweight)
}

#' @rdname setup.biteweight
#' @method setup.biteweight null
#' @param wt optional vector of biting weights, if `NULL` use 1 for all strata
#' @export
setup.biteweight.null <- function(type, model, wt = NULL, ...) {
  if (is.null(wt)) {
    wt <- rep(1, length(model$human$H))
  }
  stopifnot(is.finite(wt))
  stopifnot(wt > 0)

  biteweight$wt <- wt
  model$biteweight <- biteweight
}

#' @title Compute biting weight
#' @description a thing
#' @param biteweight an object from [MicroMoB::setup.biteweight]
#' @param t time
#' @export
compute.biteweight <- function(biteweight, t) {
  stopifnot(is.finite(t))
  UseMethod("compute.biteweight", biteweight)
}

#' @rdname compute.biteweight
#' @method compute.biteweight null
#' @export
compute.biteweight.null <- function(biteweight, t) {
  return(biteweight$wt)
}
