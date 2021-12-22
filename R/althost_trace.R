#' @title Setup trace driven alternative blood hosts
#' @description This model complies with the visitors component interface. It adds
#' a named list `model$alternative`.
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @param O a vector of length `p`, or a matrix with `p` rows and `tmax` columns,
#' or `NULL` to set to `0` (no visitors)
#' @export
setup_alternative_trace <- function(model, O = NULL) {
  stopifnot(inherits(model, "MicroMoB"))
  p <- model$global$p
  tmax <- model$global$tmax

  if (!is.null(O)) {
    if (inherits(O, "matrix")) {
      stopifnot(nrow(O) == p)
      stopifnot(ncol(O) == tmax)
    } else {
      stopifnot(length(O) == p)
      if (p > 1) {
        O <- replicate(tmax, O)
      } else {
        O <- matrix(data = O, nrow = p, ncol = tmax)
      }
    }
  } else {
    O <- matrix(data = 0, nrow = p, ncol = tmax)
  }

  model$alternative <- structure(list(), class = "trace")
  model$alternative$O <- O

}


#' @title Compute available alternative blood hosts for trace model (\eqn{O})
#' @inheritParams compute_O
#' @export
compute_O.trace <- function(model) {
  model$alternative$O[, model$global$tnow]
}
