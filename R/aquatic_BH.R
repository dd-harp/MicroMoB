# a Beverton-Holt style model of aquatic stages

#' @title Setup aquatic (immature) mosquito model with Beverton-Holt dynamics
#' @description A single compartment for all aquatic stages is modeled which
#' suffers density dependent mortality like the Beverton-Holt model.
#' @details All parameters can be passed either as a vector of length equal to `p`, a matrix with `p` rows
#' and `tmax` columns, or a matrix with `p` rows and `365` columns.
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @param stochastic should the model update deterministically or stochastically?
#' @param molt proportion of immature stages which will mature and emerge as adults each day (may be time and patch varying see [MicroMoB::time_patch_varying_parameter])
#' @param surv daily survival probability (may be time and patch varying see [MicroMoB::time_patch_varying_parameter])
#' @param K carrying capacity (may be time and patch varying see [MicroMoB::time_patch_varying_parameter])
#' @param L initial number of immature mosquitoes
#' @return no return value
#' @export
setup_aqua_BH <- function(model, stochastic, molt, surv, K, L) {
  stopifnot(inherits(model, "MicroMoB"))

  tmax <- model$global$tmax
  p <- model$global$p

  molt_mat <- time_patch_varying_parameter(param = molt, p = p, tmax = tmax)
  surv_mat <- time_patch_varying_parameter(param = surv, p = p, tmax = tmax)
  K_mat <- time_patch_varying_parameter(param = K, p = p, tmax = tmax)

  stopifnot(is.finite(molt_mat))
  stopifnot(molt_mat >= 0)
  stopifnot(is.finite(surv_mat))
  stopifnot(surv_mat >= 0)
  stopifnot(is.finite(K_mat))
  stopifnot(K_mat >= 0)

  stopifnot(length(L) == p)
  A <- rep(0, p)
  eggs <- rep(0, p)

  aqua_class <- c("BH")
  if (stochastic) {
    aqua_class <- c(aqua_class, "BH_stochastic")
    storage.mode(L) <- "integer"
    storage.mode(A) <- "integer"
    storage.mode(eggs) <- "integer"
  } else {
    aqua_class <- c(aqua_class, "BH_deterministic")
  }

  model$aqua <- structure(list(), class = aqua_class)
  model$aqua$molt <- molt_mat
  model$aqua$surv <- surv_mat
  model$aqua$K <- K_mat

  model$aqua$L <- L
  model$aqua$A <- A
  model$aqua$eggs <- eggs

}


# step function

#' @title Update aquatic (immature) mosquito populations for Beverton-Holt dynamics
#' @description This function dispatches on the second class attribute of `model$aqua`
#' for stochastic or deterministic behavior.
#' @inheritParams step_aqua
#' @return no return value
#' @export
step_aqua.BH <- function(model) {
  NextMethod()
}

#' @title Update aquatic (immature) mosquito populations for deterministic Beverton-Holt dynamics
#' @description Run a deterministic state update.
#' @inheritParams step_aqua
#' @return no return value
#' @export
step_aqua.BH_deterministic <- function(model) {

  tnow <- model$global$tnow

  molt <- model$aqua$molt[, tnow]
  surv <- model$aqua$surv[, tnow]
  K <- model$aqua$K[, tnow]
  L <- model$aqua$L

  # get eggs from adult mosquitoes
  eggs <- compute_oviposit(model)

  model$aqua$L = eggs + (1-molt)*surv*L*(K / (L + K))
  model$aqua$A = molt*surv*L*(K / (L + K))
}

#' @title Update aquatic (immature) mosquito populations for stochastic Beverton-Holt dynamics
#' @description Run a stochastic state update.
#' @inheritParams step_aqua
#' @return no return value
#' @importFrom stats rbinom
#' @export
step_aqua.BH_stochastic <- function(model) {

  tnow <- model$global$tnow
  p <- model$global$p

  molt <- model$aqua$molt[, tnow]
  surv <- model$aqua$surv[, tnow]
  K <- model$aqua$K[, tnow]
  L <- model$aqua$L

  # get eggs from adult mosquitoes
  eggs <- compute_oviposit(model)

  # survivors
  survived <- rbinom(n = p, size = L, prob = surv*(K / (L + K)))
  emerging <- rbinom(n = p, size = survived, prob = molt)

  model$aqua$L = eggs + survived - emerging
  model$aqua$A = emerging
}


# get emerging adults

#' @title Compute number of newly emerging adults from Beverton-Holt dynamics
#' @description This function dispatches on the second class attribute of `model$aqua`
#' for stochastic or deterministic behavior.
#' @inheritParams compute_emergents
#' @return a vector of length `p` giving the number of newly emerging adult in each patch
#' @export
compute_emergents.BH <- function(model) {
  model$aqua$A
}
