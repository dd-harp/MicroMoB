




#' @title Setup a human object
#' @description Setup a human object. The `model` object will have a list named `human`
#' added to it. Within `human` the `H` vector is a vector of populations, in each strata,
#' in each patch. This means that if there are 2 strata and 3 patches, the first 3
#' elements give the distribution of people in strata 1 across the 3 patches.
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

#' @rdname setup.human
#' @method setup.human strata
#' @param H a vector of human population sizes
#' @param J a matrix whose columns assign human strata to patches (rows); the
#' columns must all sum to one.
#' @export
setup.human.strata <- function(type, model, H, J = NULL, ...) {

  stopifnot(length(H) > 0)
  stopifnot(is.finite(H))
  stopifnot(H >= 0)

  if (is.null(J)) {
    J <- diag(length(H))
  }

  p <- nrow(J)
  n <- length(H)

  stopifnot(n == ncol(J))
  stopifnot(p > 0)

  if (is_binary(J)) {

    pop$J <- J
    pop$H <- H

  } else {

    stopifnot(colSums(J) == 1)

    np <- n * p

    # J does not directly assign H
    assignment_indices <- matrix(data = 1:np, nrow = p, ncol = n, byrow = FALSE)

    pop$J <- matrix(0, nrow = p, ncol = np)
    for (j in 1:p) {
      pop$J[j, assignment_indices[j, ]] <- 1
    }
    pop$H <- as.vector(J %*% diag(H))

  }

  model$human <- pop
}




#' @title Setup a time spent model
#' @description Setup a time spent model. The model object must have already
#' been initialized with a human object (see [MicroMoB::setup.human]). This adds
#' a [list] to `model` named "tisp" (Time spent).
#' @section Matrix time spent model:
#' When using `type = "matrix_day"` the model will assume that `theta` is a
#' matrix whose rows sum to \eqn{\le 1} giving the daily time spent distribution.
#' @param type a character in `c("matrix_day")`
#' @param model a model object (an [environment])
#' @param ... other arguments to be passed to type methods
#' @export
setup.timespent <- function(type, model, ...) {
  stopifnot(inherits(model, "environment"))
  stopifnot(!is.null(model$human))
  stopifnot(type %in% c("matrix_day"))
  tisp <- structure(list(), class = type)
  UseMethod("setup.timespent", tisp)
}

#' @rdname setup.timespent
#' @method setup.timespent matrix_day
#' @param theta a time spent matrix, if `NULL` the identity matrix is used
#' (everyone stays at their home patch)
#' @export
setup.timespent.matrix_day <- function(type, model, theta = NULL, ...) {

  if (is.null(theta)) {
    theta <- diag(length(model$human$H))
  } else {
    stopifnot(is.finite(theta))
    stopifnot(rowSums(theta) <= 1)
    stopifnot(rowSums(theta) > 1)
  }

  tisp$theta <- theta
  model$tisp <- tisp
}




#' @title Setup a time spent model
#' @description Setup a time spent model. The model object must have already
#' been initialized with a human object (see [MicroMoB::setup.human]). This adds
#' a [list] to `model` named "tisp" (Time spent).
#' @section Matrix time spent model:
#' When using `type = "matrix_day"` the model will assume that `theta` is a
#' matrix whose rows sum to \eqn{\le 1} giving the daily time spent distribution
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

setup.biteweight.null <- function(type, model, wt = NULL, ...) {
  if (is.null(wt)) {
    wt <- rep(1, length(model$human$H))
  }
  stopifnot(is.finite(wt))
  stopifnot(wt > 0)

  biteweight$wt <- wt
  model$biteweight <- biteweight
}

compute.biteweight <- function(biteweight, t) {
  stopifnot(is.finite(t))
  UseMethod("compute.biteweight", biteweight)
}

compute.biteweight.null <- function(biteweight, t) {
  return(biteweight$wt)
}
