# interface for aquatic (immature) mosquito populations: any model of immature mosquitoes must implement these functions

# step (update)

#' @title Update aquatic (immature) mosquito populations
#' @description This method dispatches on the type of `model$aqua`
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @return no return value
#' @export
step_aqua <- function(model) {
  UseMethod("step_aqua", model$aqua)
}

# retrieve output

#' @title Get output for aquatic (immature) mosquito populations
#' @description This method dispatches on the type of `model$aqua`. It returns
#' the current state of the aquatic component.
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @return a [data.frame]
#' @export
output_aqua <- function(model) {
  UseMethod("output_aqua", model$aqua)
}

# get emergents

#' @title Compute number of newly emerging adults (\eqn{\lambda})
#' @description This method dispatches on the type of `model$aqua`
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @return a vector of length `p` giving the number of newly emerging adult in each patch
#' @export
compute_emergents <- function(model) {
  UseMethod("compute_emergents", model$aqua)
}
