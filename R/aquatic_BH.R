# a Beverton-Holt style model of aquatic stages

#' @title Setup aquatic (immature) mosquito model with Beverton-Holt dynamics
#' @description A single compartment for all aquatic stages is modeled which
#' suffers density dependent mortality like the Beverton-Holt model.
#' @details All parameters can be passed either as a vector of length equal to `p`, a matrix with `p` rows
#' and `tmax` columns, or a matrix with `p` rows and `365` columns.
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @param molt proportion of immature stages which will mature and emerge as adults each day (may be time and patch varying see [MicroMoB::time_patch_varying_parameter])
#' @param surv daily survival probability (may be time and patch varying see [MicroMoB::time_patch_varying_parameter])
#' @param K carrying capacity (may be time and patch varying see [MicroMoB::time_patch_varying_parameter])
#' @param stochastic should the model update deterministically or stochastically?
#' @export
setup_aqua_BH <- function(model, molt, surv, K, stochastic) {
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

  aqua_class <- c("BH")
  if (stochastic) {
    aqua_class <- c(aqua_class, "BH_stochastic")
  } else {
    aqua_class <- c(aqua_class, "BH_deterministic")
  }

  model$aqua <- structure(list(), class = aqua_class)
  model$aqua$molt <- molt_mat
  model$aqua$surv <- surv_mat
  model$aqua$K <- K_mat

}
