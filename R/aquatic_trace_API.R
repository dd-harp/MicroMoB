# functions to modify parameter values of the trace (forced) aquatic model

#' @title Set daily emergence for trace (forced) aquatic mosquito model
#' @description Change the daily emergence parameter `lambda` for some times
#' and places. The parameter `lambda` is stored internally as a matrix so that `times`
#' and `places` are used to modify a submatrix, therefore the new value `lambda` should
#' either be a scalar value to update the entire submatrix or a matrix of `places` rows
#' and `times` columns.
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @param lambda new emergence
#' @param times vector of times to set the new values
#' @param places vector of places to set the new values
#' @return no return value
#' @export
set_lambda_aqua_trace <- function(model, lambda, times, places) {
  stopifnot(is.numeric(lambda), is.finite(lambda), lambda >= 0)
  if (length(lambda) == 1L) {
    lambda <- matrix(data = lambda, nrow = length(places), ncol = length(times))
  }
  stopifnot(is.matrix(lambda), nrow(lambda) == length(places), ncol(lambda) == length(times))
  stopifnot(all(times <= model$global$tmax), all(times > 0))
  stopifnot(all(places <= model$global$p), all(places > 0))

  model$aqua$lambda[places, times] <- lambda
}

#' @title Get daily emergence for Beverton-Holt aquatic mosquito model
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @param times vector of times to get values
#' @param places vector of places to get values
#' @return a [matrix]
#' @export
get_lambda_aqua_trace <- function(model, times, places) {
  return(model$aqua$lambda[places, times])
}
