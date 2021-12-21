# classes and methods to deal with immature mosquitoes and oviposition

#' @title Setup aquatic (immature) mosquito model with trace (forced) emergence
#' @description Emergence is passed as a (possibly time varying) parameter which is
#' decoupled from the adult mosquito dynamics.
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @param lambda either a vector of length equal to `p`, a matrix with `p` rows
#' and `tmax` columns, or a matrix with `p` rows and `365` columns
#' @param stochastic should the model update deterministically or stochastically?
#' @export
setup_aqua_trace <- function(model, lambda, stochastic) {
  stopifnot(inherits(model, "MicroMoB"))

  tmax <- model$global$tmax
  p <- model$global$p

  stopifnot(is.finite(lambda))
  stopifnot(lambda >= 0)

  if (inherits(lambda, "matrix")) {
    stopifnot(nrow(lambda) == p)
    if (ncol(lambda) == 365L) {
      ix <- (1:tmax) %% 365L
      ix[which(ix == 0L)] <- 365L
      lambda_mat <- lambda[, ix]
    } else if (ncol(lambda) == tmax) {
      lambda_mat <- lambda
    } else {
      stop("incorrect dimensions of lambda matrix")
    }
  } else {
    stopifnot(length(lambda) == p)
    lambda_mat <- replicate(n = tmax, expr = lambda)
  }

  aqua_class <- c("trace")
  if (stochastic) {
    aqua_class <- c(aqua_class, "trace_stochastic")
  } else {
    aqua_class <- c(aqua_class, "trace_deterministic")
  }

  model$aqua <- structure(list(), class = aqua_class)
  model$aqua$lambda <- lambda_mat

}


# step function

#' @title Update aquatic (immature) mosquito populations for forced emergence
#' @description This function does nothing as trace models are do not have
#' endogenous dynamics.
#' @inheritParams step_aqua
#' @export
step_aqua.trace <- function(model) {invisible()}


# get emerging adults

#' @title Compute number of newly emerging adults from forcing term
#' @description This function dispatches on the second argument of `model$aqua`
#' for stochastic or deterministic behavior.
#' @inheritParams compute_emergents
#' @details see [MicroMoB::compute_emergents.trace_deterministic] and [MicroMoB::compute_emergents.trace_stochastic]
#' @export
compute_emergents.trace <- function(model) {
  NextMethod()
}

#' @title Compute number of newly emerging adults from forcing term (deterministic)
#' @description Return the column of the lambda matrix for this day.
#' @inheritParams compute_emergents
#' @export
compute_emergents.trace_deterministic <- function(model) {
  return(model$aqua$lambda[, model$global$tnow])
}

#' @title Compute number of newly emerging adults from forcing term (stochastic)
#' @description Draw a Poisson distributed number of emerging adults with mean parameter
#' from the column of the trace matrix for this day.
#' @inheritParams compute_emergents
#' @importFrom stats rpois
#' @export
compute_emergents.trace_stochastic <- function(model) {
  return(rpois(n = model$global$p, lambda = model$aqua$lambda[, model$global$tnow]))
}


# add oviposition

#' @title Add eggs from oviposition to forced aquatic model
#' @description This function does nothing as trace models are not affected by
#' endogenous dynamics.
#' @inheritParams add_oviposit
#' @export
add_oviposit.trace <- function(model, eggs) {invisible()}


# compute clutch (eggs/patch/day)

#' @title Compute number of eggs laid from oviposition for each patch
#' @description This is to be used with modeling emergence as a trace.
#' @inheritParams compute_oviposit
#' @export
compute_oviposit.trace <- function(model) {invisible()}
