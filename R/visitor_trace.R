#' @title Setup trace driven visitors
#' @description This model complies with the visitors component interface. It adds
#' a named list `model$visitor`.
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @param Wd a time varying trace of visitor host availability passed to [MicroMoB::time_patch_varying_parameter]
#' or `NULL` to set to `0` (no visitors)
#' @param xd a time varying trace of visitor net infectiousness passed to [MicroMoB::time_patch_varying_parameter]
#' or `NULL` to set to `0` (no visitors)
#' @export
setup_visitor_trace <- function(model, Wd = NULL, xd = NULL) {
  stopifnot(inherits(model, "MicroMoB"))

  p <- model$global$p
  tmax <- model$global$tmax

  if (!is.null(Wd)) {
    Wd <- time_patch_varying_parameter(param = Wd, p = p, tmax = tmax)
  } else {
    Wd <- time_patch_varying_parameter(param = rep(0, p), p = p, tmax = tmax)
  }

  if (!is.null(xd)) {
    xd <- time_patch_varying_parameter(param = xd, p = p, tmax = tmax)
  } else {
    xd <- time_patch_varying_parameter(param = rep(0, p), p = p, tmax = tmax)
  }

  model$visitor <- structure(list(), class = "trace")
  model$visitor$Wd <- Wd
  model$visitor$xd <- xd

}


#' @title Compute available visitors for trace model (\eqn{W_{\delta}})
#' @inheritParams compute_Wd
#' @export
compute_Wd.trace <- function(model) {
  model$visitor$Wd[, model$global$tnow]
}

#' @title Compute net infectiousness of visitors for trace model (\eqn{x_{\delta}})
#' @inheritParams compute_xd
#' @export
compute_xd.trace <- function(model) {
  model$visitor$xd[, model$global$tnow]
}

