# observation processes for human prevalence

#' @title Update human population
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

# parameters <- list(sens = 0.95, spec = 0.91, testprop = c(0.5, 0.8, 0.9, 0.7, 0.6))
#
# n = 5
# S <- rpois(n = n, lambda = 50)
# X <- rpois(n = n, lambda = 10)
#
# H <- S + X
#
# test <- rbinom(n = n, size = H, prob = parameters$testprop)
#
# S_samp <- rhyper(nn = n, m = S, n = X, k = test)
# X_samp <- test - S_samp
#
# tp <- rbinom(n = n, size = X_samp, prob = parameters$sens)
# fn <- X_samp - tp
#
# tn <- rbinom(n = n, size = S_samp, prob = parameters$spec)
# fp <- S_samp - tn
#
# results["pos", "pos", ] <- tp
# results["neg", "neg", ] <- tn
#
# results["pos", "neg", ] <- fn
# results["neg", "pos", ] <- fp

observe_pfpr.SIS <- function(model, parameters) {

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

observe_pfpr.MOI
