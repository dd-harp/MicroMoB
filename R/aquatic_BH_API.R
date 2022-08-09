# functions to modify parameter values of the BH aquatic model

#' @title Set daily maturation probability for Beverton-Holt aquatic mosquito model
#' @description Change the daily maturation probability parameter `molt` for some times
#' and places. The parameter `molt` is stored internally as a matrix so that `times`
#' and `places` are used to modify a submatrix, therefore the new value `molt` should
#' either be a scalar value to update the entire submatrix or a matrix of `places` rows
#' and `times` columns.
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @param molt new daily maturation probability
#' @param times vector of times to set the new values
#' @param places vector of places to set the new values
#' @return no return value
#' @export
set_molt_aqua_BH <- function(model, molt, times, places) {
  stopifnot(is.numeric(molt), is.finite(molt), molt >= 0)
  if (length(molt) == 1L) {
    molt <- matrix(data = molt, nrow = length(places), ncol = length(times))
  }
  stopifnot(is.matrix(molt), nrow(molt) == length(places), ncol(molt) == length(times))
  stopifnot(all(times <= model$global$tmax), all(times > 0))
  stopifnot(all(places <= model$global$p), all(places > 0))

  model$aqua$molt[places, times] <- molt
}

#' @title Get daily maturation probability for Beverton-Holt aquatic mosquito model
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @param times vector of times to get values
#' @param places vector of places to get values
#' @return a [matrix]
#' @export
get_molt_aqua_BH <- function(model, times, places) {
  return(model$aqua$molt[places, times])
}

#' @title Set daily survival probability for Beverton-Holt aquatic mosquito model
#' @description Change the daily survival probability parameter `surv` for some times
#' and places. The parameter `surv` is stored internally as a matrix so that `times`
#' and `places` are used to modify a submatrix, therefore the new value `surv` should
#' either be a scalar value to update the entire submatrix or a matrix of `places` rows
#' and `times` columns.
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @param surv new daily survival probability
#' @param times vector of times to set the new values
#' @param places vector of places to set the new values
#' @return no return value
#' @export
set_surv_aqua_BH <- function(model, surv, times, places) {
  stopifnot(is.numeric(surv), is.finite(surv), surv >= 0)
  if (length(surv) == 1L) {
    surv <- matrix(data = surv, nrow = length(places), ncol = length(times))
  }
  stopifnot(is.matrix(surv), nrow(surv) == length(places), ncol(surv) == length(times))
  stopifnot(all(times <= model$global$tmax), all(times > 0))
  stopifnot(all(places <= model$global$p), all(places > 0))

  model$aqua$surv[places, times] <- surv
}

#' @title Get daily survival probability for Beverton-Holt aquatic mosquito model
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @param times vector of times to get values
#' @param places vector of places to get values
#' @return a [matrix]
#' @export
get_surv_aqua_BH <- function(model, times, places) {
  return(model$aqua$surv[places, times])
}

#' @title Set carrying capacity for Beverton-Holt aquatic mosquito model
#' @description Change the carrying capacity parameter `K` for some times
#' and places. The parameter `K` is stored internally as a matrix so that `times`
#' and `places` are used to modify a submatrix, therefore the new value `K` should
#' either be a scalar value to update the entire submatrix or a matrix of `places` rows
#' and `times` columns.
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @param K new carrying capacity
#' @param times vector of times to set the new values
#' @param places vector of places to set the new values
#' @return no return value
#' @export
set_K_aqua_BH <- function(model, K, times, places) {
  stopifnot(is.numeric(K), is.finite(K), K >= 0)
  if (length(K) == 1L) {
    K <- matrix(data = K, nrow = length(places), ncol = length(times))
  }
  stopifnot(is.matrix(K), nrow(K) == length(places), ncol(K) == length(times))
  stopifnot(all(times <= model$global$tmax), all(times > 0))
  stopifnot(all(places <= model$global$p), all(places > 0))

  model$aqua$K[places, times] <- K
}

#' @title Get carrying capacity for Beverton-Holt aquatic mosquito model
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @param times vector of times to get values
#' @param places vector of places to get values
#' @return a [matrix]
#' @export
get_K_aqua_BH <- function(model, times, places) {
  return(model$aqua$K[places, times])
}
