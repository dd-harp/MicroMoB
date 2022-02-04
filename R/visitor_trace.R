#' @title Setup trace driven visitors
#' @description This model complies with the visitors component interface. It adds
#' a named list `model$visitor`.
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @param Wd a time varying trace of visitor host availability passed to [MicroMoB::time_patch_varying_parameter]
#' or `NULL` to set to `0` (no visitors)
#' @param xd a time varying trace of visitor net infectiousness passed to [MicroMoB::time_patch_varying_parameter]
#' or `NULL` to set to `0` (no visitors)
#' @return no return value
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


#' @title Get parameters for trace driven visitors
#' @description The JSON config file should have two entries:
#'  * Wd: vector or matrix (see [MicroMoB::time_patch_varying_parameter] for valid dimensions)
#'  * xd: vector or matrix (see [MicroMoB::time_patch_varying_parameter] for valid dimensions)
#'
#' For interpretation of the entries, please read [MicroMoB::setup_visitor_trace].
#' @param path a file path to a JSON file
#' @return a named [list]
#' @importFrom jsonlite read_json
#' @examples
#' # to see an example of proper JSON input, run the following
#' library(jsonlite)
#' par <- list(
#'  "Wd" = rep(1, 5),
#'  "xd" = rep(0.01, 365)
#' )
#' toJSON(par)
#' @export
get_config_visitor_trace <- function(path) {
  pars <- read_json(path = file.path(path), simplifyVector = TRUE)

  stopifnot(length(pars) == 2L)

  stopifnot(is.numeric(pars$Wd))
  stopifnot(is.vector(pars$Wd) | is.matrix(pars$Wd))

  stopifnot(is.numeric(pars$xd))
  stopifnot(is.vector(pars$xd) | is.matrix(pars$xd))

  return(pars)
}


#' @title Compute available visitors for trace model (\eqn{W_{\delta}})
#' @inheritParams compute_Wd
#' @return a vector of length `p` giving biting availability of visitors at each patch
#' @export
compute_Wd.trace <- function(model) {
  model$visitor$Wd[, model$global$tnow]
}

#' @title Compute net infectiousness of visitors for trace model (\eqn{x_{\delta}})
#' @inheritParams compute_xd
#' @return a vector of length `p` giving net infectiousness of visitors at each patch
#' @export
compute_xd.trace <- function(model) {
  model$visitor$xd[, model$global$tnow]
}

