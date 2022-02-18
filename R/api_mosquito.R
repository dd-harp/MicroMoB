# functions for setting up APIs, includes global methods and those for aquatic/adult components (see inst/plumber/mosquito/plumber.R for API)

#' @noRd
#' @importFrom jsonlite read_json
api_config_global <- function(path) {
  pars <- read_json(path = file.path(path), simplifyVector = TRUE)

  # mosy-only
  if (length(pars) == 6L) {
    stopifnot(is.numeric(pars$tmax))
    stopifnot(is.numeric(pars$p))

    stopifnot(is.character(pars$aqua_path))
    stopifnot(file.exists(pars$aqua_path))
    stopifnot(pars$aqua_model %in% c("BH", "trace"))

    stopifnot(is.character(pars$adult_path))
    stopifnot(file.exists(pars$adult_path))
    stopifnot(pars$adult_model %in% c("RM"))
    # full sim
  } else if (length(pars) == 10L) {
    stop("configuration of full transmission simulation has not been implemented yet")
  } else {
    stop("invalid global config file provided")
  }

  return(pars)
}


#' @noRd
#' @importFrom jsonlite unbox
api_setup_global_parameters <- function(path, res) {

  global_pars <- api_config_global(path = path)

  parenv <- parent.frame()

  parenv$parameters <- list()
  parenv$parameters$global <- global_pars

  res$status <- 200
  return(list(msg = unbox("global parameters successfully read in")))
}



#' @noRd
#' @importFrom jsonlite unbox
api_setup_aqua_parameters <- function(res) {

  parenv <- parent.frame()

  if (!exists(x = "parameters", envir = parenv)) {
    res$status <- 400
    return(list(error = unbox("model parameters not yet specified")))
  }

  global_pars <- parenv$parameters$global

  if (global_pars$aqua_model == "BH") {
    aqua_pars <- get_config_aqua_BH(path = global_pars$aqua_path)
  } else if (global_pars$aqua_model == "trace") {
    aqua_pars <- get_config_aqua_trace(path = global_pars$aqua_path)
  } else {
    res$status <- 400
    return(list(error = unbox("invalid aquatic model specified")))
  }

  parenv$parameters$aqua <- aqua_pars

  res$status <- 200
  return(list(msg = unbox("aquatic component parameters successfully read in")))
}


#' @noRd
#' @importFrom jsonlite unbox
api_setup_adult_parameters <- function(res) {

  parenv <- parent.frame()

  if (!exists(x = "parameters", envir = parenv)) {
    res$status <- 400
    return(list(error = unbox("model parameters not yet specified")))
  }

  global_pars <- parenv$parameters$global

  if (global_pars$adult_model == "RM") {
    adult_pars <- get_config_mosquito_RM(path = global_pars$adult_path)
  } else {
    res$status <- 400
    return(list(error = unbox("invalid adult model specified")))
  }

  parenv$parameters$adult <- adult_pars

  res$status <- 200
  return(list(msg = unbox("adult component parameters successfully read in")))
}


#' @noRd
#' @importFrom jsonlite unbox
api_get_parameters_global <- function(res) {
  if (!exists(x = "parameters", where = parent.frame())) {
    res$status <- 400
    return(list(error = unbox("model parameters not yet specified")))
  } else {
    parenv <- parent.frame()
    return(parenv$parameters$global)
  }
}


#' @noRd
#' @importFrom jsonlite unbox
api_get_parameters_adult <- function(res) {
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
api_get_parameters_aqua <- function(res) {
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
api_setup_model_object <- function(res) {

  parenv <- parent.frame()

  if (!exists(x = "parameters", envir = parenv)) {
    res$status <- 400
    return(list(error = unbox("model parameters not yet specified")))
  }

  global_pars <- parenv$parameters$global

  # setup global model obj
  parenv$mod <- make_MicroMoB(tmax = global_pars$tmax, p = global_pars$p)

  res$status <- 200
  return(list(msg = unbox("model object successfully set up")))
}


#' @noRd
#' @importFrom jsonlite unbox
api_setup_aqua <- function(res) {

  parenv <- parent.frame()

  if (!exists(x = "parameters", envir = parenv)) {
    res$status <- 400
    return(list(error = unbox("model parameters not yet specified")))
  }

  global_pars <- parenv$parameters$global
  aqua_pars <-  parenv$parameters$aqua

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
  return(list(msg = unbox("aquatic component successfully set up")))
}


#' @noRd
#' @importFrom jsonlite unbox
api_setup_adult <- function(res) {

  parenv <- parent.frame()

  if (!exists(x = "parameters", envir = parenv)) {
    res$status <- 400
    return(list(error = unbox("model parameters not yet specified")))
  }

  global_pars <- parenv$parameters$global
  adult_pars <-  parenv$parameters$adult

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

  res$status <- 200
  return(list(msg = unbox("adult component successfully set up")))
}




#' @noRd
#' @importFrom jsonlite unbox
api_step_aqua <- function(res) {

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
api_step_adult <- function(res) {

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
api_get_output_aqua <- function(res) {

  parenv <- parent.frame()

  if (!exists(x = "mod", envir = parenv)) {
    res$status <- 400
    return(list(error = unbox("model object not found")))
  }

  output_aqua(parenv$mod)
}


#' @noRd
api_get_output_adult <- function(res) {

  parenv <- parent.frame()

  if (!exists(x = "mod", envir = parenv)) {
    res$status <- 400
    return(list(error = unbox("model object not found")))
  }

  output_mosquitoes(parenv$mod)
}

#' @noRd
api_clock_tick <- function(res) {

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
