# classes and methods to deal with immature mosquitoes and oviposition

#' @title Setup aquatic (immature) mosquito model with trace (forced) emergence
#' @description Emergence is passed as a (possibly time varying) parameter which is
#' decoupled from the adult mosquito dynamics. This module assumes `l` and `p`
#' are equivalent, as emergence rates are given for `p`.
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @param lambda daily emergence of mosquitoes, may be time and patch varying, see [MicroMoB::time_patch_varying_parameter]
#' @param stochastic should the model update deterministically or stochastically?
#' @return no return value
#' @export
setup_aqua_trace <- function(model, lambda, stochastic) {
  stopifnot(inherits(model, "MicroMoB"))

  tmax <- model$global$tmax
  l <- model$global$l

  stopifnot(is.finite(lambda))
  stopifnot(lambda >= 0)

  lambda_mat <- time_patch_varying_parameter(param = lambda, p = l, tmax = tmax)

  aqua_class <- c("trace")
  if (stochastic) {
    aqua_class <- c(aqua_class, "trace_stochastic")
  } else {
    aqua_class <- c(aqua_class, "trace_deterministic")
  }

  model$aqua <- structure(list(), class = aqua_class)
  model$aqua$lambda <- lambda_mat

}


#' @title Get parameters for aquatic (immature) model with forced emergence
#' @description The JSON config file should have two entries:
#'  * stochastic: a boolean value
#'  * lambda: a scalar, vector, or matrix (row major). It will be passed to
#'  [MicroMoB::time_patch_varying_parameter], see that function's documentation for
#'  appropriate dimensions.
#'
#' For interpretation of the entries, please read [MicroMoB::setup_aqua_trace].
#' @param path a file path to a JSON file
#' @return a named [list]
#' @importFrom jsonlite read_json
#' @examples
#' # to see an example of proper JSON input, run the following
#' library(jsonlite)
#' t <- 10 # number of days to simulate
#' par <- list(
#'  "stochastic" = FALSE,
#'  "lambda" = rpois(n = t, lambda = 10)
#' )
#' toJSON(par, pretty = TRUE)
#' @export
get_config_aqua_trace <- function(path) {
  pars <- read_json(path = file.path(path), simplifyVector = TRUE)
  stopifnot(length(pars) == 2L)
  stopifnot(is.logical(pars$stochastic))
  stopifnot(length(pars$stochastic) == 1L)
  stopifnot(is.numeric(pars$lambda))
  stopifnot(is.vector(pars$lambda) | is.matrix(pars$lambda))
  return(pars)
}


# output

#' @title Get output for aquatic (immature) mosquito populations with forced emergence
#' @description This function returns an empty [data.frame] as trace models do
#' not have endogenous dynamics.
#' @inheritParams output_aqua
#' @return a [data.frame]
#' @export
output_aqua.trace <- function(model) {data.frame()}


# step function

#' @title Update aquatic (immature) mosquito populations for forced emergence
#' @description This function does nothing as trace models do not have
#' endogenous dynamics.
#' @inheritParams step_aqua
#' @return no return value
#' @export
step_aqua.trace <- function(model) {invisible()}


# get emerging adults

#' @title Compute number of newly emerging adults from forcing term
#' @description This function dispatches on the second class attribute of `model$aqua`
#' for stochastic or deterministic behavior.
#' @inheritParams compute_emergents
#' @details see [MicroMoB::compute_emergents.trace_deterministic] and [MicroMoB::compute_emergents.trace_stochastic]
#' @return no return value
#' @export
compute_emergents.trace <- function(model) {
  NextMethod()
}

#' @title Compute number of newly emerging adults from forcing term (deterministic)
#' @description Return the column of the lambda matrix for this day.
#' @inheritParams compute_emergents
#' @return a vector of length `l` giving the number of newly emerging adult in each patch
#' @export
compute_emergents.trace_deterministic <- function(model) {
  return(model$aqua$lambda[, model$global$tnow])
}

#' @title Compute number of newly emerging adults from forcing term (stochastic)
#' @description Draw a Poisson distributed number of emerging adults with mean parameter
#' from the column of the trace matrix for this day.
#' @inheritParams compute_emergents
#' @return a vector of length `l` giving the number of newly emerging adult in each patch
#' @importFrom stats rpois
#' @export
compute_emergents.trace_stochastic <- function(model) {
  return(rpois(n = model$global$l, lambda = model$aqua$lambda[, model$global$tnow]))
}
