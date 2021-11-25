# pathogen models

#' @title Setup a pathogen model
#' @description Setup a time spent model. The model object must have already
#' been initialized with a human object (see [MicroMoB::setup_human]). The `model` object will have a list named `pathogen`
#' added to it.
#' @seealso [MicroMoB::setup_pathogen.rm]
#' @param type a character in `c("rm")`
#' @param model a model object (an [environment])
#' @param ... other arguments to be passed to type methods
#' @export
setup_pathogen <- function(type, model, ...) {
  stopifnot(inherits(model, "environment"))
  stopifnot(!is.null(model$human))
  pathogen <- structure(list(), class = type)
  UseMethod("setup_pathogen", pathogen)
}

#' @title Setup daily time spent model
#' @description A model of where people spend their time over the day.
#' @inheritParams setup_pathogen
#' @param rm_parameters a list of parameters
#' @export
setup_pathogen.rm <- function(type, model, rm_parameters, ...) {

  stopifnot(length(rm_parameters[["r"]]) == 1L)
  stopifnot(length(rm_parameters[["b"]]) == 1L)
  stopifnot(length(rm_parameters[["c"]]) == 1L)

  stopifnot(is.finite(rm_parameters[["r"]]))
  stopifnot(is.finite(rm_parameters[["b"]]))
  stopifnot(is.finite(rm_parameters[["c"]]))

  n <- model$human$n
  X <- rm_parameters[["X"]]

  stopifnot(length(X) == n)
  stopifnot(X >= 0)
  stopifnot(X <= 1)

  pathogen$r <- rm_parameters[["r"]]
  pathogen$b <- rm_parameters[["b"]]
  pathogen$c <- rm_parameters[["c"]]
  pathogen$X <- X

  model$pathogen <- pathogen
}


#' @export
setup_pathogen.default <- function(type, model, ...) {
  stop("setup_pathogen has no method for dispatch type ", type)
}


# get prevalence

#' @title Compute resident net infectiousness (\eqn{x})
#' @description \eqn{x}, the net infectiousness of residents, is the probability
#' that a mosquito would become infected after biting a resident It is the prevalence
#' multiplied by the transmission efficiency.
#' @param pathogen an object from [MicroMoB::setup_pathogen]
#' @param t time
#' @seealso [MicroMoB::compute_x.rm]
#' @export
compute_x <- function(pathogen, t) {
  UseMethod("compute_x", pathogen)
}

#' @title Compute simple visitor availability (\eqn{x})
#' @inheritParams compute_x
#' @return a vector of dimension \eqn{p \times 1}
#' @export
compute_x.rm <- function(pathogen, t) {
  return(pathogen$X * pathogen$c)
}

#' @export
compute_x.default <- function(pathogen, t) {
  stop("compute_x has no method for dispatch type ", class(pathogen))
}
