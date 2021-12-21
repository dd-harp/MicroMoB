# interface for humans: any model of humans must implement these functions

# update (step)

#' @title Update human population
#' @description This method dispatches on the type of `model$human`.
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @export
step_humans <- function(model) {
  UseMethod("step_humans", model$human)
}


#' @title Compute available humans
#' @description This is normally computed as \deqn{W = \Psi^{\intercal} \cdot w_{f} H}.
#' This method dispatches on the type of `model$human`.
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @return a vector of length `p`
#' @export
compute_W <- function(model) {
  UseMethod("compute_W", model$human)
}


#' @title Compute net infectiousness of humans
#' @description In a Ross-Macdonald style transmission model, this is computed as
#' \deqn{x = c X} This method dispatches on the type of `model$human`.
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @return a vector of length `n`
#' @export
compute_x <- function(model) {
  UseMethod("compute_x", model$human)
}
