#' @title Setup trace driven visitors
#' @description This model complies with the visitors component interface. It adds
#' a named list `model$visitor`.
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @param Wd a vector of length `p`, or a matrix with `p` rows and `tmax` columns,
#' or `NULL` to set to `0` (no visitors)
#' @param xd a vector of length `p`, or a matrix with `p` rows and `tmax` columns,
#' or `NULL` to set to `0` (no visitors)
#' @export
setup_visitor_trace <- function(model, Wd = NULL, xd = NULL) {
  stopifnot(inherits(model, "MicroMoB"))
  p <- model$global$p
  tmax <- model$global$tmax

  if (!is.null(Wd)) {
    if (inherits(Wd, "matrix")) {
      stopifnot(nrow(Wd) == p)
      stopifnot(ncol(Wd) == tmax)
    } else {
      stopifnot(length(Wd) == p)
      if (p > 1) {
        Wd <- replicate(tmax, Wd)
      } else {
        Wd <- matrix(data = Wd, nrow = p, ncol = tmax)
      }
    }
  } else {
    Wd <- matrix(data = 0, nrow = p, ncol = tmax)
  }

  if (!is.null(xd)) {
    if (inherits(xd, "matrix")) {
      stopifnot(nrow(xd) == p)
      stopifnot(ncol(xd) == tmax)
    } else {
      stopifnot(length(xd) == p)
      if (p > 1) {
        xd <- replicate(tmax, xd)
      } else {
        xd <- matrix(data = xd, nrow = p, ncol = tmax)
      }
    }
  } else {
    xd <- matrix(data = 0, nrow = p, ncol = tmax)
  }

  model$visitor <- structure(list(), class = "trace")
  model$visitor$Wd <- Wd
  model$visitor$xd <- xd

}


#' @title Compute available visitors for trace model (\eqn{W_{\delta}})
#' @inheritParams compute_Wd
#' @export
compute_Wd.trace <- function(model) {
  model$visitor$Wd[, model$global$tnow]
}

#' @title Compute net infectiousness of visitors for trace model (\eqn{x_{\delta}})
#' @inheritParams compute_xd
#' @export
compute_xd.trace <- function(model) {
  model$visitor$xd[, model$global$tnow]
}

