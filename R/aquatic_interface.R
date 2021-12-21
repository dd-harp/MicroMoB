# interface for aquatic (immature) mosquito populations: any model of immature mosquitoes must implement these functions

# step (update)

#' @title Update aquatic (immature) mosquito populations
#' @description This method dispatches on the type of `model$aqua`
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @export
step_aqua <- function(model) {
  UseMethod("step_aqua", model$aqua)
}

# get emergents

#' @title Compute number of newly emerging adults
#' @description This method dispatches on the type of `model$aqua`
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @export
compute_emergents <- function(model) {
  UseMethod("compute_emergents", model$aqua)
}


# add oviposition

#' @title Add eggs from oviposition to aquatic model
#' @description This method dispatches on the type of `model$aqua`
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @param eggs a vector of length `p` giving eggs for each place
#' @export
add_oviposit <- function(model, eggs) {
  UseMethod("add_oviposit", model$aqua)
}


# compute oviposition (eggs/patch/day)

#' @title Compute number of eggs laid from oviposition for each patch
#' @description This method dispatches on the type of `model$aqua`
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @export
compute_oviposit <- function(model) {
  UseMethod("compute_oviposit", model$aqua)
}
