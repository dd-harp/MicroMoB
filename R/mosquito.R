#' @title Setup a mosquito object
#' @description Setup a mosquito object. The `model` object will have a list named `mosquito`
#' added to it.
#' @seealso [MicroMoB::setup_mosquito.rm_day]
#' @param type a character in `c("rm_day")`
#' @param model a model object (an [environment])
#' @param ... other arguments to be passed to type methods
#' @export
setup_mosquito <- function(type, model, ...) {
  stopifnot(inherits(model, "environment"))
  pop <- structure(list(), class = type)
  UseMethod("setup_human", pop)
}

#' @title Setup a Ross-Macdonald mosquito object
#' @description a model
#' @inheritParams setup_mosquito
#' @param mosquito_rm_pars parameters for Ross-Macdonald mosquito model
#' @export
setup_mosquito.rm_day <- function(type, model, mosquito_rm_pars, ...) {

}

#' @export
setup_mosquito.default <- function(type, model, ...) {
  stop("setup_human has no method for dispatch type ", type)
}
