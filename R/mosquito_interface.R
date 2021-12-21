# interface for mosquitoes: any model of mosquitoes must implement these functions

#' @title Update mosquito population
#' @description This method dispatches on the type of `model$mosquito`
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @export
step_mosquitoes <- function(model) {
  UseMethod("step_mosquitoes", model$mosquito)
}

#' @title Compute mosquito feeding rate (\eqn{f})
#' @description This method dispatches on the type of `model$mosquito`
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @param B a vector of length `p` giving total blood host availability by patch
#' @return a vector of length `p`
#' @export
compute_f <- function(model, B) {
  UseMethod("compute_f", model$mosquito)
}

#' #' @title Compute human blood feeding fraction (\eqn{q})
#' #' @description This method dispatches on the type of `model$mosquito`
#' #' @param model an object from [MicroMoB::make_MicroMoB]
#' #' @return a vector of length `p`
#' #' @export
#' compute_q <- function(model) {
#'   UseMethod("compute_q", model$mosquito)
#' }

#' @title Compute density of infective mosquitoes (\eqn{Z})
#' @description This method dispatches on the type of `model$mosquito`
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @return a vector of length `p`
#' @export
compute_Z <- function(model) {
  UseMethod("compute_Z", model$mosquito)
}
