# observation processes for human prevalence

#' @title Observe PfPR in human strata
#' @description This method dispatches on the type of `model$human`.
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @param parameters a named [list], should have elements `sens` (sensitivity),
#' `spec` (specificity), and a vector of length equal to number of strata `testprop`
#' which gives the proportion of each strata to be tested.
#' @return an [array] of counts, with actual condition as first dimension and tested condition
#' as the second dimension, and the third dimension is the human strata
#' @export
observe_pfpr <- function(model, parameters) {
  UseMethod("observe_pfpr", model$human)
}

#' @title Observe PfPR in human strata for SIS model
#' @inheritParams observe_pfpr
#' @return an [array] of counts, with actual condition as first dimension and tested condition
#' as the second dimension, and the third dimension is the human strata
#' @importFrom stats rbinom rhyper
#' @export
observe_pfpr.SIS <- function(model, parameters) {

  n <- model$global$n

  stopifnot(length(parameters$testprop) == n)

  results <- array(data = 0, dim = c(2, 2, n), dimnames = list(c("pos", "neg"), c("pos", "neg"), NULL))

  test <- rbinom(n = n, size = model$human$H, prob = parameters$testprop)

  # if sampling randomly, how many true S and I are sampled for testing?
  S <- model$human$H - model$human$X
  S_samp <- rhyper(nn = n, m = S, n = model$human$X, k = test)
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


#' @title Observe PfPR in human strata for SIP model
#' @inheritParams observe_pfpr
#' @return an [array] of counts, with actual condition as first dimension and tested condition
#' as the second dimension, and the third dimension is the human strata
#' @importFrom stats rbinom rhyper
#' @export
observe_pfpr.SIP <- function(model, parameters) {

  n <- model$global$n

  stopifnot(length(parameters$testprop) == n)

  results <- array(data = 0, dim = c(2, 2, n), dimnames = list(c("pos", "neg"), c("pos", "neg"), NULL))

  H <- rowSums(model$human$SIP)
  test <- rbinom(n = n, size = H, prob = parameters$testprop)

  # how many of these tests are randomly allocated to infected persons?
  I_samp <- rhyper(nn = n, m = model$human$SIP[, "I"], n = rowSums(model$human$SIP[, c("S", "P")]), k = test)
  SP_samp <- test - I_samp

  tp <- rbinom(n = n, size = I_samp, prob = parameters$sens)
  fn <- I_samp - tp

  tn <- rbinom(n = n, size = SP_samp, prob = parameters$spec)
  fp <- SP_samp - tn

  results["pos", "pos", ] <- tp
  results["neg", "neg", ] <- tn

  results["pos", "neg", ] <- fn
  results["neg", "pos", ] <- fp

  return(results)
}
