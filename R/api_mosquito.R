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

hello_world <- function() {
  "Hello, world!"
}

#' @title Run Plumber API for mosquito-only simulation
#' @param ... arguments passed to [plumber::pr_run]
#' @importFrom plumber plumb_api pr_run
#' @export
run_api_mosquito <- function(...) {
  pr_run(pr = plumb_api(package = "MicroMoB", name = "mosquito", edit = FALSE), ...)
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


#' @noRd
#' @importFrom jsonlite unbox
put_model_object_mosquito <- function(res) {

  parenv <- parent.frame()

  if (!exists(x = "parameters", envir = parenv)) {
    res$status <- 400
    return(list(error = unbox("model parameters not yet specified")))
  }

  global_pars <- parenv$parameters$global
  adult_pars <- parenv$parameters$adult
  aqua_pars <- parenv$parameters$aqua

  # setup global model obj
  parenv$mod <- make_MicroMoB(tmax = global_pars$tmax, p = global_pars$p)

  # setup adult component
  if (global_pars$adult_model == "RM") {
    setup_mosquito_RM(
      model = parenv$mod,
      stochastic = adult_pars$stochastic,
      f = adult_pars$f,
      q = adult_pars$q,
      eip = adult_pars$eip,
      p = adult_pars$p,
      psi = adult_pars$psi,
      nu = adult_pars$nu,
      M = adult_pars$M,
      Y = adult_pars$Y,
      Z = adult_pars$Z
    )
  } else {
    res$status <- 400
    return(list(error = unbox("unknown adult mosquito model specified")))
  }

  # setup aquatic component
  if (global_pars$aqua_model == "trace") {
    setup_aqua_trace(
      model = parenv$mod,
      lambda = aqua_pars$lambda,
      stochastic = aqua_pars$stochastic
    )
  } else if (global_pars$aqua_model == "BH") {
    setup_aqua_BH(
      model = parenv$mod,
      stochastic = aqua_pars$stochastic,
      molt = aqua_pars$molt,
      surv = aqua_pars$surv,
      K = aqua_pars$K,
      L = aqua_pars$L
    )
  } else {
    res$status <- 400
    return(list(error = unbox("unknown aquatic mosquito model specified")))
  }

  res$status <- 200
  return(list(msg = unbox("model successfully set up")))
}


#' @noRd
put_step_aqua <- function(res) {

  parenv <- parent.frame()

  if (!exists(x = "mod", envir = parenv)) {
    res$status <- 400
    return(list(error = unbox("model object not found")))
  }

  if (parenv$mod$global$tnow > parenv$mod$global$tmax) {
    res$status <- 400
    return(list(error = unbox("attempt to step model out of defined time span")))
  }

  step_aqua(parenv$mod)
}


#' @noRd
put_step_adult <- function(res) {

  parenv <- parent.frame()

  if (!exists(x = "mod", envir = parenv)) {
    res$status <- 400
    return(list(error = unbox("model object not found")))
  }

  if (parenv$mod$global$tnow > parenv$mod$global$tmax) {
    res$status <- 400
    return(list(error = unbox("attempt to step model out of defined time span")))
  }

  step_mosquitoes(parenv$mod)
}


#' @noRd
get_output_aqua <- function(res) {

  parenv <- parent.frame()

  if (!exists(x = "mod", envir = parenv)) {
    res$status <- 400
    return(list(error = unbox("model object not found")))
  }

  output_aqua(parenv$mod)
}


#' @noRd
get_output_adult <- function(res) {

  parenv <- parent.frame()

  if (!exists(x = "mod", envir = parenv)) {
    res$status <- 400
    return(list(error = unbox("model object not found")))
  }

  output_mosquitoes(parenv$mod)
}

#' @noRd
clock_tick <- function(res) {

  parenv <- parent.frame()

  if (!exists(x = "mod", envir = parenv)) {
    res$status <- 400
    return(list(error = unbox("model object not found")))
  }

  if (parenv$mod$global$tnow > parenv$mod$global$tmax) {
    res$status <- 400
    return(list(error = unbox("attempt to step model out of defined time span")))
  }

  parenv$mod$global$tnow <- parenv$mod$global$tnow + 1L

}
