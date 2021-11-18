


#' @title Setup other blood hosts
#' @description Setup a other blood hosts (vertebrate mammals) object.
#' The `model` object will have a list named `otherhosts` added to it.
#' @section strata: see [MicroMoB::setup.human.strata]
#' @param type a character in `c("simple")`
#' @param model a model object (an [environment])
#' @param ... other arguments to be passed to type methods
#' @export
setup.otherhosts <- function(type, model, ...) {
  stopifnot(inherits(model, "environment"))
  stopifnot(!is.null(model$human))
  stopifnot(type %in% c("simple"))
  x <- structure(list(), class = type)
  UseMethod("setup.otherhosts", x)
}

setup.otherhosts.simple <- function(type, model, B = NULL, ...) {
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


compute_B.otherhosts <- function(otherhosts, t) {
  UseMethod("compute_B.otherhosts", otherhosts)
}

compute_B.otherhosts.simple <- function(otherhosts, t) {
  return(otherhosts$B)
}
