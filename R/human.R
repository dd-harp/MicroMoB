# setup human model object

#' @title Setup a human object
#' @description Setup a human object. The `model` object will have a list named `human`
#' added to it.
#' @seealso [MicroMoB::setup_human.strata]
#' @param type a character in `c("strata")`
#' @param model a model object (an [environment])
#' @param ... other arguments to be passed to type methods
#' @export
setup_human <- function(type, model, ...) {
  stopifnot(inherits(model, "environment"))
  pop <- structure(list(), class = type)
  UseMethod("setup_human", pop)
}

#' @title Setup a human model with strata
#' @description This sets up a human model object. If you need help getting
#' lumped data into the correct format, see [MicroMoB::strata_to_residency_proportion].
#' @inheritParams setup_human
#' @param H a vector of human population sizes
#' @param J a matrix whose columns assign human strata to patches (rows); the
#' columns must all sum to one. If `J` is `NULL` then a diagonal matrix of ones
#' is used; this assumes that the only level of population stratification is
#' patch of residence.
#' @export
setup_human.strata <- function(type, model, H, J = NULL, ...) {

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

#' @export
setup_human.default <- function(type, model, ...) {
  stop("setup_human has no method for dispatch type ", type)
}


# setup biting weights

#' @title Setup a biting weight model
#' @description Setup a biting weight model. The model object must have already
#' been initialized with a human object (see [MicroMoB::setup_human]). This adds
#' a [list] to `model$human` named `biteweight`.
#' @seealso [MicroMoB::setup_biteweight.simple]
#' @param type a character in `c("simple")`
#' @param model a model object (an [environment])
#' @param ... other arguments to be passed to type methods
#' @export
setup_biteweight <- function(type, model, ...) {
  stopifnot(inherits(model, "environment"))
  stopifnot(!is.null(model$human))
  biteweight <- structure(list(), class = type)
  UseMethod("setup_biteweight", biteweight)
}

#' @title Setup simple biting weight model
#' @description A simple vector of biting weights that do not change over time.
#' @inheritParams setup_biteweight
#' @param wf optional vector of biting weights, if `NULL` use 1 for all strata
#' @export
setup_biteweight.simple <- function(type, model, wf = NULL, ...) {
  if (is.null(wf)) {
    wf <- rep(1, length(model$human$n))
  }
  stopifnot(is.finite(wf))
  stopifnot(wf > 0)

  biteweight$wf <- wf
  model$human$biteweight <- biteweight
}

#' @export
setup_biteweight.default <- function(type, model, ...) {
  stop("setup_biteweight has no method for dispatch type ", type)
}


#' @title Compute biting weight
#' @description Return the biting weights for humans at a given time.
#' @seealso [MicroMoB::compute_wf.simple]
#' @param biteweight an object from [MicroMoB::setup_biteweight]
#' @param t time
#' @export
compute_wf <- function(biteweight, t) {
  stopifnot(is.finite(t))
  UseMethod("compute_wf", biteweight)
}

#' @title Compute biting weight for simple model
#' @inheritParams compute_wf
#' @export
compute_wf.simple <- function(biteweight, t) {
  return(biteweight$wf)
}

#' @export
compute_wf.default <- function(biteweight, t) {
  stop("compute_wf has no method for dispatch type ", class(biteweight))
}
