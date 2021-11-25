#' @title Setup a mosquito object
#' @description Setup a mosquito object. The `model` object will have a list named `mosquito`
#' added to it.
#' @seealso [MicroMoB::setup_mosquito.trace]
#' @param type a character in `c("rm", "null)`
#' @param model a model object (an [environment])
#' @param ... other arguments to be passed to type methods
#' @export
setup_mosquito <- function(type, model, ...) {
  stopifnot(inherits(model, "environment"))
  pop <- structure(list(), class = type)
  UseMethod("setup_mosquito", pop)
}

#' @title Setup a trace mosquito object
#' @description This model simply returns a predefined vector `Z`, giving the number
#' of infectious female mosquitoes ine each patch during each day.
#' @inheritParams setup_mosquito
#' @param Z a matrix whose columns are days and rows are patches
#' @param xi a
#' @export
setup_mosquito.trace <- function(type, model, Z, xi, ...) {
  stopifnot(inherits(Z, "matrix"))
  stopifnot(is.finite(Z))
  stopifnot(Z >= 0)

  stopifnot(length(xi) > 0)

  if (length(xi) == 1L) {
    pop$xi <- 1
  } else {
    stopifnot(is.finite(xi))
    stopifnot(sum(xi) == 1)
    pop$xi <- xi
  }
  pop$Z <- Z

  model$mosquito <- pop
}

#' @export
setup_mosquito.default <- function(type, model, ...) {
  stop("setup_mosquito has no method for dispatch type ", type)
}


# get number of infectious mosquitoes

#' @title Compute infectious mosquitoes (\eqn{Z})
#' @param mosquito an object from [MicroMoB::setup_mosquito]
#' @param t time
#' @seealso [MicroMoB::compute_Z.trace]
#' @export
compute_Z <- function(mosquito, t) {
  UseMethod("compute_Z", mosquito)
}

#' @title Compute simple visitor availability (\eqn{W_{\delta}})
#' @inheritParams compute_Z
#' @return a vector of dimension \eqn{p \times 1}
#' @export
compute_Z.trace <- function(mosquito, t) {
  return(mosquito$Z[, t, drop = FALSE])
}

#' @export
compute_Z.default <- function(mosquito, t) {
  stop("compute_Z has no method for dispatch type ", class(mosquito))
}
