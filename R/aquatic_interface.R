# interface for aquatic (immature) mosquito populations: any model of immature mosquitoes must implement these functions

# step (update)

#' @title Update aquatic (immature) mosquito populations
#' @description This method dispatches on the type of `model$aqua`
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @export
step_aqua <- function(model) {
  UseMethod("step_aqua", model$aqua)
}

# get emergents

#' @title Compute number of newly emerging adults (\eqn{\lambda})
#' @description This method dispatches on the type of `model$aqua`
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @export
compute_emergents <- function(model) {
  UseMethod("compute_emergents", model$aqua)
}
