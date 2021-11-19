is_binary <- function(x) {
  uniq <- unique(as.vector(x))
  if (all(length(uniq) == 2L)) {
    return(all(uniq %in% c(0, 1)))
  } else {
    return(FALSE)
  }
}

#' @title Functional form type 0
#' @description An equation of the form
#' \deqn{f(x) = \frac{f_{x} s_{f} x}{1 + s_{f} x}}
#' which has an asymptote at \eqn{f_{x}} and \eqn{s_{f}} controls the rate at which
#' the asymptote is reached.
#' @param x the value
#' @param fx the asymptote
#' @param sf controls slope
Ff_0 <- function(x, fx = 1, sf = 1){
  stopifnot(is.finite(x))
  fx * (sf * x / (1 + sf * x))
}

#' @title Helper function for lumped population strata (proportional assignment)
#' @description If input is given as a vector of population sizes per-strata, lumped
#' over patches, and a separate matrix whose columns describe how each strata is
#' distributed over patches, this function calculates the residency matrix and
#' population size for the overall stratification of both residency and strata.
#' @param H_strata a vector of population size by strata
#' @param J_strata a matrix chose columns sum to one giving the distribution of
#' strata populations over patches
#' @return a [list] with three elements:
#'  * `assignment_indices`: provides a mapping from patch (rows) and strata (columns)
#'     into the "unrolled" vector `H`
#'  * `J`: the residency matrix mapping elements in `H` to patches
#'  * `H`: the overall population distribution over strata and patches
#' @examples
#' # taken from package tests
#' J <- matrix(
#'    c(0.3, 0.5, 0.2,
#'    0.1, 0.6, 0.3), nrow = 3, ncol = 2, byrow = FALSE
#' )
#' H <- c(50, 60)
#' H_overall <- J %*% diag(H)
#' residency <- strata_to_residency_proportion(H_strata = H, J_strata = J)
#' @export
strata_to_residency_proportion <- function(H_strata, J_strata) {

  stopifnot(inherits(J_strata, "matrix"))
  stopifnot(is.finite(H_strata))
  stopifnot(H_strata >= 0)
  stopifnot(length(H_strata) == ncol(J_strata))
  stopifnot(colSums(J_strata) == 1)

  p <- nrow(J_strata)
  s <- ncol(J_strata)
  n <- p*s

  # we need to determine the actual residency matrix
  pop <- list()

  pop$J <- do.call(what = cbind, args = replicate(n = s, expr = diag(p), simplify = FALSE))
  pop$H <- as.vector(J_strata %*% diag(H_strata))

  return(pop)
}

#' @title Helper function for lumped population strata (counts)
#' @description If input is given as a matrix of population counts per strata (columns)
#' and patch (rows), this function calculates the residency matrix and
#' population size for the overall stratification of both residency and strata.
#' @param H_counts a matrix of population counts
#' @return a [list] with three elements:
#'  * `J`: the residency matrix mapping elements in `H` to patches
#'  * `H`: the overall population distribution over strata and patches
#' @examples
#' # taken from package tests
#' J <- matrix(
#'    c(0.3, 0.5, 0.2,
#'    0.1, 0.6, 0.3), nrow = 3, ncol = 2, byrow = FALSE
#' )
#' H <- c(50, 60)
#' H_overall <- J %*% diag(H)
#' residency <- strata_to_residency_proportion(H_strata = H, J_strata = J)
#' @export
strata_to_residency_counts <- function(H_counts) {

  stopifnot(inherits(H_counts, "matrix"))
  stopifnot(is.finite(H_counts))
  stopifnot(H_counts >= 0)

  p <- nrow(H_counts)
  s <- ncol(H_counts)

  pop <- list()

  pop$H <- as.vector(H_counts)
  pop$J <- do.call(what = cbind, args = replicate(n = s, expr = diag(p), simplify = FALSE))

  return(pop)
}
