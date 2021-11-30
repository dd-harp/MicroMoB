# pathogen models

#' @title Setup a pathogen model
#' @description Setup a time spent model. The model object must have already
#' been initialized with a human object (see [MicroMoB::setup_human]). The `model` object will have a list named `pathogen`
#' added to it.
#' @seealso [MicroMoB::setup_pathogen.rm]
#' @param type a character in `c("rm")`
#' @param model a model object (from [MicroMoB::setup_model_object])
#' @param ... other arguments to be passed to type methods
#' @export
setup_pathogen <- function(type, model, ...) {
  stopifnot(inherits(model, "micro_mob"))
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

  stopifnot(length(rm_parameters[["r"]]) == 1L | length(rm_parameters[["r"]]) == model$gloabl$n)
  stopifnot(length(rm_parameters[["b"]]) == 1L | length(rm_parameters[["b"]]) == model$gloabl$n)
  stopifnot(length(rm_parameters[["c"]]) == 1L | length(rm_parameters[["c"]]) == model$gloabl$n)

  stopifnot(is.logical(rm_parameters[["stochastic"]]))
  stopifnot(is.finite(rm_parameters[["r"]]))
  stopifnot(is.finite(rm_parameters[["b"]]))
  stopifnot(is.finite(rm_parameters[["c"]]))

  n <- model$human$n
  X <- rm_parameters[["X"]]

  stopifnot(length(X) == n)
  stopifnot(X >= 0) # X is density or count

  pathogen$r <- rm_parameters[["r"]]
  pathogen$b <- rm_parameters[["b"]]
  pathogen$c <- rm_parameters[["c"]]
  pathogen$stochastic <- rm_parameters[["stochastic"]]

  if (length(pathogen$b) == 1L) {
    pathogen$b <- diag(rep(pathogen$b, model$gloabl$n))
  } else {
    pathogen$b <- diag(pathogen$b)
  }

  if (length(pathogen$r) == 1L) {
    pathogen$r <- rep(pathogen$r, model$gloabl$n)
  }

  if (length(pathogen$c) == 1L) {
    pathogen$c <- rep(pathogen$c, model$gloabl$n)
  }

  pathogen$X <- X
  pathogen$h <- rep(0, n)

  model$pathogen <- pathogen
}


#' @export
setup_pathogen.default <- function(type, model, ...) {
  stop("setup_pathogen has no method for dispatch type ", type)
}


# get force of infection

#' @title Compute force of infection on human hosts (\eqn{h})
#' @description \eqn{h}, the force of infection on each human strata, is calculated
#' after calling [MicroMoB::compute_biting] which computes the biting distribution matrix
#' and feeding rates.
#' @param model a model object (an object from [MicroMoB::setup_model_object])
#' @param t time
#' @export
compute_h <- function(model, t) {

  Z <- compute_Z(model$mosquito)
  b <- model$pathogen$b
  beta <- model$biting$beta
  f <- model$biting$f
  q <- model$biting$q

  model$pathogen$h <- (b %*% beta) %*% (f*q*Z)

}


# get prevalence

#' @title Compute resident net infectiousness (\eqn{x})
#' @description \eqn{x}, the net infectiousness of residents, is the probability
#' that a mosquito would become infected after biting a resident It is the prevalence
#' multiplied by the transmission efficiency.
#' @param model an object from [MicroMoB::setup_model_object]
#' @param t time
#' @seealso [MicroMoB::compute_x.rm]
#' @export
compute_x <- function(model, t) {
  UseMethod("compute_x", model$pathogen)
}

#' @title Compute resident net infectiousness for Ross-Macdonland pathogen model (\eqn{x})
#' @inheritParams compute_x
#' @return a vector of dimension \eqn{p \times 1}
#' @export
compute_x.rm <- function(model, t) {
  X <- model$pathogen$X / model$human$H
  stopifnot(is.finite(X))
  return(X * model$pathogen$c)
}

#' @export
compute_x.default <- function(model, t) {
  stop("compute_x has no method for dispatch type ", class(model$pathogen))
}


# time step

#' @title Update the pathogen model
#' @param model a model object (an object from [MicroMoB::setup_model_object])
#' @param t time
#' @seealso [MicroMoB::step_pathogen.rm]
#' @export
step_pathogen <- function(model, t) {
  UseMethod("step_pathogen", model$pathogen)
}

#' @title Update the Ross-Macdonald pathogen model
#' @inheritParams step_pathogen
#' @importFrom stats pexp rbinom
#' @export
step_pathogen.rm <- function(model, t) {
  S <- model$human$H - model$pathogen$X
  X <- model$pathogen$X

  if (model$pathogen$stochastic) {
    new_infections <- rbinom(n = length(S), size = S, prob = pexp(q = model$pathogen$h))
    recoveries <- rbinom(n = length(X), size = X, prob = pexp(q = model$pathogen$r))
  } else {
    new_infections <- S*pexp(q = model$pathogen$h)
    recoveries <- X*pexp(q = model$pathogen$r)
  }

  model$pathogen$X <- model$pathogen$X + new_infections - recoveries
  stopifnot(model$pathogen$X >= 0)
}

#' @export
step_pathogen.default <- function(model, t) {
  stop("step_pathogen has no method for dispatch type ", class(model$pathogen))
}

