# plumber API for mosquito-only simulation

#* Setup parameters for mosquito-only simulation
#* @put /config_mosquito
put_config_mosquito <- function(path, res) {

  global_pars <- MicroMoB::get_config_mosquito_MicroMoB(path = path)

  if (global_pars$aqua_model == "BH") {
    aqua_pars <- MicroMoB::get_config_aqua_BH(path = global_pars$aqua_path)
  } else if (global_pars$aqua_model == "trace") {
    aqua_pars <- MicroMoB::get_config_aqua_trace(path = global_pars$aqua_path)
  } else {
    res$status <- 400
    return(list(error = jsonlite::unbox("invalid aquatic model specified")))
  }

  if (global_pars$adult_model == "RM") {
    adult_pars <- MicroMoB::get_config_mosquito_RM(path = global_pars$adult_path)
  } else {
    res$status <- 400
    return(list(error = jsonlite::unbox("invalid adult model specified")))
  }

  assign(x = "parameters", value = list(), pos = .GlobalEnv)
  .GlobalEnv$parameters$global <- global_pars
  .GlobalEnv$parameters$aqua <- aqua_pars
  .GlobalEnv$parameters$adult <- adult_pars

  res$status <- 200
  return(list(msg = jsonlite::unbox("model parameters successfully read in")))
}

#* Get parameters for adult mosquitoes for mosquito-only simulation
#* @get /parameters_adults
#* @serializer json list(pretty = TRUE)
get_parameters_adult_mosquito <- function(res) {
  if (!exists(x = "parameters", where = .GlobalEnv)) {
    res$status <- 400
    return(list(error = jsonlite::unbox("model parameters not yet specified")))
  } else {
    return(.GlobalEnv$parameters$adult)
  }
}


#* Get parameters for aquatic (immature) mosquitoes for mosquito-only simulation
#* @get /parameters_adults
#* @serializer json list(pretty = TRUE)
get_parameters_aqua_mosquito <- function(res) {
  if (!exists(x = "parameters", where = .GlobalEnv)) {
    res$status <- 400
    return(list(error = jsonlite::unbox("model parameters not yet specified")))
  } else {
    return(.GlobalEnv$parameters$aqua)
  }
}
