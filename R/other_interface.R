# interface for visitors: any model of visitors must implement these functions

#' @title Compute available alternative blood hosts (\eqn{O})
#' @description This method dispatches on the type of `model$alternative`.
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @return a vector of length `p`
#' @export
compute_O <- function(model) {
  UseMethod("compute_O", model$alternative)
}
