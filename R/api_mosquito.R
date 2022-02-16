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


#' @title Run Plumber API for mosquito-only simulation
#' @param ... arguments passed to [plumber::pr_run]
#' @importFrom plumber plumb_api pr_run
#' @export
run_api_mosquito <- function(...) {
  pr_run(pr = plumb_api(package = "MicroMoB", name = "mosquito", edit = FALSE), ...)
}
