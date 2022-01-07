# interface for visitors: any model of visitors must implement these functions

#' @title Compute available visitors (\eqn{W_{\delta}})
#' @description This method dispatches on the type of `model$visitor`.
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @return a vector of length `p` giving biting availability of visitors at each patch
#' @export
compute_Wd <- function(model) {
  UseMethod("compute_Wd", model$visitor)
}

#' @title Compute net infectiousness of visitors (\eqn{x_{\delta}})
#' @description This method dispatches on the type of `model$visitor`.
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @return a vector of length `p` giving net infectiousness of visitors at each patch
#' @export
compute_xd <- function(model) {
  UseMethod("compute_xd", model$visitor)
}
