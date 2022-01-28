# interface for mosquitoes: any model of mosquitoes must implement these functions

#' @title Update mosquito population
#' @description This method dispatches on the type of `model$mosquito`
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @return no return value
#' @export
step_mosquitoes <- function(model) {
  UseMethod("step_mosquitoes", model$mosquito)
}

#' @title Compute mosquito feeding rate (\eqn{f})
#' @description This method dispatches on the type of `model$mosquito`
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @param B a vector of length `p` giving total blood host availability by patch
#' @return a vector of length `p` giving the per-capita blood feeding rate of mosquitoes in each patch
#' @export
compute_f <- function(model, B) {
  UseMethod("compute_f", model$mosquito)
}

#' @title Compute human blood feeding fraction (\eqn{q})
#' @description This method dispatches on the type of `model$mosquito`
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @param W a vector of length `p` giving human availability by patch (\eqn{W})
#' @param Wd a vector of length `p` giving visitor availability by patch (\eqn{W_{\delta}})
#' @param B a vector of length `p` giving total blood host availability by patch (\eqn{B})
#' @return a vector of length `p` giving the proportion of bites taken on human hosts in each patch
#' @export
compute_q <- function(model, W, Wd, B) {
  UseMethod("compute_q", model$mosquito)
}

#' @title Compute density of infective mosquitoes (\eqn{Z})
#' @description This method dispatches on the type of `model$mosquito`. \eqn{Z}
#' is also known as the "sporozoite rate" in malariaology.
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @return a vector of length `p` giving the density of infected and infectious mosquitoes in each patch
#' @export
compute_Z <- function(model) {
  UseMethod("compute_Z", model$mosquito)
}


# compute oviposition (eggs/patch/day)

#' @title Compute number of eggs laid from oviposition for each patch
#' @description This method dispatches on the type of `model$mosquito`
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @return a vector of length `p` giving the total number of eggs laid by adult mosquitoes in each patch
#' @export
compute_oviposit <- function(model) {
  UseMethod("compute_oviposit", model$mosquito)
}
