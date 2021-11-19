# utility functions

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
#' lumped data into the correct format, see [MicroMoB::strata_to_residency].
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
#' a [list] to `model$human` named `biteweight` (Time spent).
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
#' @param wt optional vector of biting weights, if `NULL` use 1 for all strata
#' @export
setup_biteweight.simple <- function(type, model, wt = NULL, ...) {
  if (is.null(wt)) {
    wt <- rep(1, length(model$human$n))
  }
  stopifnot(is.finite(wt))
  stopifnot(wt > 0)

  biteweight$wt <- wt
  model$human$biteweight <- biteweight
}

#' @export
setup_biteweight.default <- function(type, model, ...) {
  stop("setup_biteweight has no method for dispatch type ", type)
}


#' @title Compute biting weight
#' @description Return the biting weights for humans at a given time
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
  return(biteweight$wt)
}

#' @export
compute_wf.default <- function(biteweight, t) {
  stop("compute_wf has no method for dispatch type ", class(biteweight))
}
