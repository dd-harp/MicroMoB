# observation processes for human prevalence

#' @title Observe prevalence in human strata
#' @description This method dispatches on the type of `model$human`.
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @param parameters a named [list], should have elements `sens` (sensitivity),
#' `spec` (specificity), and a vector of length equal to number of strata `testprop`
#' which gives the proportion of each strata to be tested.
#' @return an [array] of counts, with actual condition as first dimension and tested condition
#' as the second dimension, and the third dimension is the human strata
#' @export
observe_prev <- function(model, parameters) {
  UseMethod("observe_prev", model$human)
}

#' @title Observe prevalence in human strata for SIS model
#' @inheritParams observe_prev
#' @return an [array] of counts, with actual condition as first dimension and tested condition
#' as the second dimension, and the third dimension is the human strata
#' @importFrom stats rbinom rhyper
#' @export
observe_prev.SIS <- function(model, parameters) {

  n <- model$global$n

  stopifnot(length(parameters$testprop) == n)

  results <- array(data = 0, dim = c(2, 2, n), dimnames = list(c("pos", "neg"), c("pos", "neg"), NULL))

  test <- rbinom(n = n, size = model$human$H, prob = parameters$testprop)

  S_samp <- rhyper(nn = n, m = model$human$S, n = model$human$X, k = test)
  X_samp <- test - S_samp

  tp <- rbinom(n = n, size = X_samp, prob = parameters$sens)
  fn <- X_samp - tp

  tn <- rbinom(n = n, size = S_samp, prob = parameters$spec)
  fp <- S_samp - tn

  results["pos", "pos", ] <- tp
  results["neg", "neg", ] <- tn

  results["pos", "neg", ] <- fn
  results["neg", "pos", ] <- fp

  return(results)
}


#' @title Observe prevalence in human strata for MOI model
#' @inheritParams observe_prev
#' @return an [array] of counts, with actual condition as first dimension and tested condition
#' as the second dimension, and the third dimension is the human strata
#' @export
observe_prev.MOI  <- function(model, parameters) {
  stop("not yet implemented for the MOI model")
}


#' @title Observe prevalence in human strata for SIR model
#' @inheritParams observe_prev
#' @return an [array] of counts, with actual condition as first dimension and tested condition
#' as the second dimension, and the third dimension is the human strata
#' @export
observe_prev.SIR  <- function(model, parameters) {
  stop("not yet implemented for the SIR model")
}


#' @title Observe prevalence in human strata for SIR model
#' @inheritParams observe_prev
#' @return an [array] of counts, with actual condition as first dimension and tested condition
#' as the second dimension, and the third dimension is the human strata
#' @export
observe_prev.SIP  <- function(model, parameters) {
  stop("not yet implemented for the SIR model")
}
