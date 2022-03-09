# this is a null model of mosquito dynamics that is only for testing/verifying aquatic models.
# it implements a single method `compute_oviposit` and all other methods throw an error

#' @title Setup null mosquito model
#' @description This is a null model of mosquito dynamics that is only for testing/verifying aquatic models.
#' It implements a single method [MicroMoB::compute_oviposit.trace] and all other methods throw an error.
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @param oviposit a vector of length `p` used as a return value for [MicroMoB::compute_oviposit]
#' @return no return value
#' @export
setup_mosquito_trace <- function(model, oviposit) {
  stopifnot(inherits(model, "MicroMoB"))
  stopifnot(length(oviposit) == model$global$p)
  stopifnot(is.finite(oviposit))
  stopifnot(oviposit >= 0)

  model$mosquito <- structure(list(), class = "trace")
  model$mosquito$oviposit <- oviposit
}

#' @title Get parameters for null mosquito model
#' @description The JSON config file should have 1 entry:
#'  * oviposit: vector
#'
#' For interpretation of the entries, please read [MicroMoB::setup_mosquito_trace].
#' @param path a file path to a JSON file
#' @return a named [list]
#' @importFrom jsonlite read_json
#' @examples
#' # to see an example of proper JSON input, run the following
#' library(jsonlite)
#' par <- list(
#'  "oviposit" = rep(1, 5)
#' )
#' toJSON(par, pretty = TRUE)
#' @export
get_config_mosquito_trace <- function(path) {
  pars <- read_json(path = file.path(path), simplifyVector = TRUE)

  stopifnot(length(pars) == 1L)

  stopifnot(is.numeric(pars$oviposit))
  stopifnot(is.vector(pars$oviposit))

  return(pars)
}


# output

#' @title Get output for null mosquito populations
#' @description This function returns an empty [data.frame] as trace models do
#' not have endogenous dynamics.
#' @inheritParams output_mosquitoes
#' @return a [data.frame]
#' @export
output_mosquitoes.trace <- function(model) {data.frame()}

#' @title Update null mosquito population
#' @inheritParams step_mosquitoes
#' @return no return value
#' @export
step_mosquitoes.trace <- function(model) {
  stop("trace adult model does not implement a step method")
}

#' @title Compute null mosquito feeding rate (\eqn{f})
#' @inheritParams compute_f
#' @return no return value
#' @export
compute_f.trace <- function(model, B) {
  stop("trace adult model does not implement this method")
}

#' @title Compute null human blood feeding fraction (\eqn{q})
#' @inheritParams compute_q
#' @return no return value
#' @export
compute_q.trace <- function(model, W, Wd, B) {
  stop("trace adult model does not implement this method")
}

#' @title Compute null density of infective mosquitoes (\eqn{Z})
#' @inheritParams compute_Z
#' @return no return value
#' @export
compute_Z.trace <- function(model) {
  stop("trace adult model does not implement this method")
}


# compute oviposition (eggs/patch/day)

#' @title Compute number of eggs laid from oviposition for each patch for null model
#' @description This method dispatches on the type of `model$mosquito`
#' @inheritParams compute_oviposit
#' @return a vector of length `p` giving the total number of eggs laid by adult mosquitoes in each patch
#' @export
compute_oviposit.trace <- function(model) {
  model$mosquito$oviposit
}
