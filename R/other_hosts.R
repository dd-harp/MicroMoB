# other vertebrate blood hosts

#' @title Setup other blood hosts
#' @description Setup a other blood hosts (vertebrate mammals) object.
#' The `model` object will have a list named `otherhosts` added to it.
#' @seealso [setup_otherhosts.simple]
#' @param type a character in `c("simple")`
#' @param model a model object (an [environment])
#' @param ... other arguments to be passed to type methods
#' @export
setup_otherhosts <- function(type, model, ...) {
  stopifnot(inherits(model, "environment"))
  stopifnot(!is.null(model$human))
  x <- structure(list(), class = type)
  UseMethod("setup_otherhosts", x)
}

#' @title Setup simple other blood hosts
#' @inheritParams setup_otherhosts
#' @param B an optional vector of search weights by patch, otherwise will
#' be set to all 0.
#' @param zeta search weights will be weighted as \eqn{B^{\Zeta}}
#' @export
setup_otherhosts.simple <- function(type, model, B = NULL, zeta = 1, ...) {
  p <- model$human$p
  if (is.null(B)) {
    x$B <- rep(0, p)
  } else {
    stopifnot(length(B) == p)
    stopifnot(is.finite(B))
    stopifnot(B >= 0)
    x$B <- B
  }

  stopifnot(zeta >= 0)
  stopifnot(length(zeta) == 1 | length(zeta) == p)
  x$zeta <- zeta

  model$otherhosts <- x
}

#' @export
setup_otherhosts.default <- function(type, model, ...) {
  stop("setup_otherhosts has no method for dispatch type ", type)
}


# compute other host availability

#' @title Compute other host availability (\eqn{B})
#' @param otherhosts an object from [MicroMoB::setup_otherhosts]
#' @param t time
#' @seealso [MicroMoB::compute_B.simple]
#' @export
compute_B <- function(otherhosts, t) {
  UseMethod("compute_B", otherhosts)
}

#' @title Compute simple other host availability (\eqn{B})
#' @inheritParams compute_B
#' @return a vector of dimension \eqn{p \times 1}
#' @export
compute_B.simple <- function(otherhosts, t) {
  return(otherhosts$B^otherhosts$zeta)
}

#' @export
compute_B.default <- function(otherhosts, t) {
  stop("compute_B has no method for dispatch type ", class(otherhosts))
}


# human visitors

#' @title Setup visitors
#' @description Setup a human visitors (vertebrate mammals) object.
#' These are persons who are present in the spatial domain but are not
#' residents of any patch (transient visitors).
#' The `model` object will have a list named `visitors` added to it.
#' @seealso [setup_visitors.simple]
#' @param type a character in `c("simple")`
#' @param model a model object (an [environment])
#' @param ... other arguments to be passed to type methods
#' @export
setup_visitors <- function(type, model, ...) {
  stopifnot(inherits(model, "environment"))
  stopifnot(!is.null(model$human))
  x <- structure(list(), class = type)
  UseMethod("setup_visitors", x)
}

#' @title Setup simple visitors
#' @inheritParams setup_visitors
#' @description If both parameters are `NULL` then
#' @param W_delta a vector of blood feeding attempt weights by patch
#' @param x_delta a vector of prevalence in visitors by patch
#' @export
setup_visitors.simple <- function(type, model, W_delta = NULL, x_delta = NULL, ...) {
  p <- model$human$p

  if (is.null(W_delta)) {
    x$W_delta <- rep(0, p)
  } else {
    stopifnot(length(W_delta) == p)
    stopifnot(is.finite(W_delta))
    stopifnot(W_delta >= 0)
    x$W_delta <- W_delta
  }

  if (is.null(x_delta)) {
    x$x_delta <- rep(0, p)
  } else {
    stopifnot(length(x_delta) == p)
    stopifnot(is.finite(x_delta))
    stopifnot(x_delta >= 0)
    x$x_delta <- x_delta
  }

  model$visitors <- x
}

#' @export
setup_visitors.default <- function(type, model, ...) {
  stop("setup_visitors has no method for dispatch type ", type)
}

#' @title Compute visitor availability (\eqn{W_{\delta}})
#' @param visitors an object from [MicroMoB::setup_visitors]
#' @param t time
#' @seealso [MicroMoB::compute_W_delta.simple]
#' @export
compute_W_delta <- function(visitors, t) {
  UseMethod("compute_W_delta", visitors)
}

#' @title Compute simple visitor availability (\eqn{W_{\delta}})
#' @inheritParams compute_W_delta
#' @return a vector of dimension \eqn{p \times 1}
#' @export
compute_W_delta.simple <- function(visitors, t) {
  return(visitors$W_delta)
}

#' @export
compute_W_delta.default <- function(visitors, t) {
  stop("compute_W_delta has no method for dispatch type ", class(visitors))
}


#' @title Compute visitor net infectiousness (\eqn{x_{\delta}})
#' @description \eqn{x_{\delta}}, the net infectiousness of visitors, is the probability
#' that a mosquito would become infected after biting a visitor. It is the prevalence
#' multiplied by the transmission efficiency.
#' @param visitors an object from [MicroMoB::setup_visitors]
#' @param t time.
#' @seealso [MicroMoB::compute_x_delta.simple]
#' @export
compute_x_delta <- function(visitors, t) {
  UseMethod("compute_x_delta", visitors)
}

#' @title Compute simple visitor net infectiousness (\eqn{x_{\delta}})
#' @inheritParams compute_x_delta
#' @return a vector of dimension \eqn{p \times 1}
#' @export
compute_x_delta.simple <- function(visitors, t) {
  return(visitors$x_delta)
}

#' @export
compute_x_delta.default <- function(visitors, t) {
  stop("compute_x_delta has no method for dispatch type ", class(visitors))
}
