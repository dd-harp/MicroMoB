#' @title Make a model object
#' @description The model object is a hashed [environment]. By default it contains
#' a single list, `model$global` storing global state.
#' @param tmax number of days to simulate
#' @param p number of places
#' @export
make_MicroMoB <- function(tmax, p) {
  stopifnot(is.finite(tmax))
  stopifnot(tmax > 0)
  stopifnot(is.finite(p))
  stopifnot(p > 0)
  object <- structure(new.env(hash = TRUE), class = "MicroMoB")
  object$global <- list(tmax = as.integer(tmax), tnow = 1L, p = as.integer(p))
  return(object)
}


