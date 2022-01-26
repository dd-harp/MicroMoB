#' @title Setup trace driven alternative blood hosts
#' @description This model complies with the visitors component interface. It adds
#' a named list `model$alternative`.
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @param O a time varying trace passed to [MicroMoB::time_patch_varying_parameter]
#' or `NULL` to set to `0` (no alternative blood hosts)
#' @return no return value
#' @export
setup_alternative_trace <- function(model, O = NULL) {
  stopifnot(inherits(model, "MicroMoB"))

  p <- model$global$p
  tmax <- model$global$tmax

  if (!is.null(O)) {
    O <- time_patch_varying_parameter(param = O, p = p, tmax = tmax)
  } else {
    O <- time_patch_varying_parameter(param = rep(0, p), p = p, tmax = tmax)
  }

  model$alternative <- structure(list(), class = "trace")
  model$alternative$O <- O

}


#' @title Compute available alternative blood hosts for trace model (\eqn{O})
#' @inheritParams compute_O
#' @return a vector of length `p` giving biting availability of other blood hosts at each patch
#' @export
compute_O.trace <- function(model) {
  model$alternative$O[, model$global$tnow]
}
