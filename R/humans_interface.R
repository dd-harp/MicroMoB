# interface for humans: any model of humans must implement these functions

# update (step)

#' @title Update human population
#' @description This method dispatches on the type of `model$human`.
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @export
step_humans <- function(model) {
  UseMethod("step_humans", model$human)
}


#' @title Compute available humans (\eqn{W})
#' @description This is normally computed as \deqn{W = \Psi^{\top} \cdot w_{f} H}.
#' This method dispatches on the type of `model$human`.
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @return a vector of length `p` giving the biting availability of human hosts at each patch
#' @export
compute_W <- function(model) {
  UseMethod("compute_W", model$human)
}


#' @title Compute human biting weights (\eqn{w_{f}})
#' @description This method dispatches on the type of `model$human`.
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @return a vector of length `n` giving the biting weights of human hosts in each stratum
#' @export
compute_wf <- function(model) {
  UseMethod("compute_wf", model$human)
}


#' @title Compute net infectiousness of humans (\eqn{x})
#' @description In a Ross-Macdonald style transmission model, this is computed as
#' \deqn{x = c X} This method dispatches on the type of `model$human`.
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @return a vector of length `n` giving the net infectiousness of human hosts in each stratum
#' @export
compute_x <- function(model) {
  UseMethod("compute_x", model$human)
}


#' @title Compute human population strata sizes (\eqn{H})
#' @description This method dispatches on the type of `model$human`.
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @return a vector of length `n` giving the size of each human population stratum
#' @export
compute_H <- function(model) {
  UseMethod("compute_H", model$human)
}


#' @title Compute time at risk matrix (\eqn{\Psi})
#' @description The time at risk matrix is \eqn{\Psi = \Theta \xi} This method dispatches on the type of `model$human`.
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @return a matrix with `n` rows and `p` columns, the time at risk matrix
#' @export
compute_Psi <- function(model) {
  UseMethod("compute_Psi", model$human)
}
