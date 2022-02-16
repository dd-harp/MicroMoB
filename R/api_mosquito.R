# helper functions for mosquito-only simulation (see inst/plumber/mosquito/plumber.R for API)

#' @title Get parameters for configuring a mosquito only Micro-MoB simulation
#' @description The JSON config file should have:
#'  * tmax: integer
#'  * p: integer
#'  * aqua_path: file path
#'  * aqua_model: character in `BH`, `trace`
#'  * adult_path: file path
#'  * adult_model: character in `RM`
#'
#' @param path a file path to a JSON file
#' @return a named [list]
#' @importFrom jsonlite read_json
#' @export
get_config_mosquito_MicroMoB <- function(path) {
  pars <- read_json(path = file.path(path), simplifyVector = TRUE)

  stopifnot(length(pars) == 6L)

  stopifnot(is.numeric(pars$tmax))
  stopifnot(is.numeric(pars$p))

  stopifnot(is.character(pars$aqua_path))
  stopifnot(file.exists(pars$aqua_path))
  stopifnot(pars$aqua_model %in% c("BH", "trace"))

  stopifnot(is.character(pars$adult_path))
  stopifnot(file.exists(pars$adult_path))
  stopifnot(pars$adult_model %in% c("RM"))

  return(pars)
}

#' @noRd
#' @importFrom jsonlite unbox
put_config_mosquito <- function(path, res) {

  global_pars <- get_config_mosquito_MicroMoB(path = path)

  if (global_pars$aqua_model == "BH") {
    aqua_pars <- get_config_aqua_BH(path = global_pars$aqua_path)
  } else if (global_pars$aqua_model == "trace") {
    aqua_pars <- get_config_aqua_trace(path = global_pars$aqua_path)
  } else {
    res$status <- 400
    return(list(error = unbox("invalid aquatic model specified")))
  }

  if (global_pars$adult_model == "RM") {
    adult_pars <- get_config_mosquito_RM(path = global_pars$adult_path)
  } else {
    res$status <- 400
    return(list(error = unbox("invalid adult model specified")))
  }

  parenv <- parent.frame()

  parenv$parameters <- list()
  parenv$parameters$global <- global_pars
  parenv$parameters$aqua <- aqua_pars
  parenv$parameters$adult <- adult_pars

  res$status <- 200
  return(list(msg = unbox("model parameters successfully read in")))
}

#' @noRd
#' @importFrom jsonlite unbox
get_parameters_adult_mosquito <- function(res) {
  if (!exists(x = "parameters", where = parent.frame())) {
    res$status <- 400
    return(list(error = unbox("model parameters not yet specified")))
  } else {
    parenv <- parent.frame()
    return(parenv$parameters$adult)
  }
}

#' @noRd
#' @importFrom jsonlite unbox
get_parameters_aqua_mosquito <- function(res) {
  if (!exists(x = "parameters", where = parent.frame())) {
    res$status <- 400
    return(list(error = unbox("model parameters not yet specified")))
  } else {
    parenv <- parent.frame()
    return(parenv$parameters$aqua)
  }
}

#' @title Run Plumber API for mosquito-only simulation
#' @param ... arguments passed to [plumber::pr_run]
#' @importFrom plumber plumb_api pr_run
#' @export
run_api_mosquito <- function(...) {
  pr_run(pr = plumb_api(package = "MicroMoB", name = "mosquito", edit = FALSE), ...)
}
