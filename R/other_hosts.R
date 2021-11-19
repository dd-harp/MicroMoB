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
#' @export
setup_otherhosts.simple <- function(type, model, B = NULL, ...) {
  p <- model$human$p
  if (is.null(B)) {
    x$B <- rep(0, p)
  } else {
    stopifnot(length(B) == p)
    stopifnot(is.finite(B))
    stopifnot(B >= 0)
    x$B <- B
  }
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
  return(otherhosts$B)
}

#' @export
compute_B.default <- function(otherhosts, t) {
  stop("compute_B has no method for dispatch type ", class(otherhosts))
}
