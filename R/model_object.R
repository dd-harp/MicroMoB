#' @title Setup a Micro-MoB model object
#'
#' @export
setup_model_object <- function() {
  model <- structure(new.env(), class = "micro_mob")
  model$biting <- structure(list(), class = "biting")
  model$global <- structure(list(), class = "global")
  return(model)
}
