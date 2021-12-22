#' @noRd
is_binary <- function(x) {
  uniq <- unique(as.vector(x))
  if (all(length(uniq) == 2L)) {
    return(all(uniq %in% c(0, 1)))
  } else {
    return(FALSE)
  }
}


#' @noRd
approx_equal <- function(a, b, tol = sqrt(.Machine$double.eps)) {
  abs(a - b) < tol
}

#' @noRd
divmod <- function(a,b){
  a <- as.integer(a)
  b <- as.integer(b)
  c(
    quo = a %/% b,
    rem = a %% b
  )
}

#' @noRd
distribute <- function(n,p){
  n <- as.integer(n)
  p <- as.integer(p)
  distn <- rep(0L,p)
  div <- divmod(n,p)
  for(i in 0:(p-1)){
    distn[i+1] <- div[["quo"]] + (i < div[["rem"]])
  }
  return(distn)
}


#' @title Sample a stochastic vector
#' @description Given a vector of counts in cells, `x` and a stochastic matrix `prob`, each
#' row of which describes a probability distribution of how that cell should be
#' distributed among bins, sample destination bins for each cell count, and return
#' a vector giving the number of counts in bins. It is conceptually similar to
#' "stochastically" distributing the vector as `x %*% prob`, which gives the
#' expectation.
#' @param x a vector
#' @param prob a matrix, it must have number of rows equal to `x` and rows that sum
#' to one
#' @return a vector of length equal to the number of columns of `prob`
#' @importFrom stats rmultinom
#' @export
sample_stochastic_vector <- function(x, prob) {
  stopifnot(length(x) == nrow(prob))
  if (ncol(prob) == 1L) {
    return(x)
  }
  samp <- vapply(X = 1:length(x), FUN = function(i) {
    rmultinom(n = 1, size = x[i], prob = prob[i, ])
  }, FUN.VALUE = numeric(ncol(prob)), USE.NAMES = FALSE)
  rowSums(samp)
}


#' @title Sample a stochastic matrix
#' @description `x` is a matrix with arbitrary number of rows but whose columns
#' are equal to the number of bins that the stochastic matrix `prob` parameterizes
#' a distribution over. Each row of `x` gives a distribution of counts over bins
#' and is resampled according to `prob`. It is conceptually similar to
#' "stochastically" distributing the matrix as `x %*% prob`, which gives the
#' expectation.
#' @param x a matrix
#' @param prob a matrix, it must have number of columns equal to the number of columns of `x` and rows that sum
#' to one
#' @return a matrix whose dimensions equal the original `x`
#' @export
sample_stochastic_matrix <- function(x, prob) {
  stopifnot(ncol(x) == ncol(prob))
  samp <- lapply(X = 1:nrow(x), FUN = function(i) {
    sample_stochastic_vector(x = x[i, ], prob = prob)
  })
  do.call(rbind, samp)
}


#' @title Helper function for lumped population strata (proportional assignment)
#' @description If input is given as a vector of population sizes per-strata, lumped
#' over patches, and a separate matrix whose columns describe how each strata is
#' distributed over patches, this function calculates the residency matrix and
#' population size for the overall stratification of both residency and strata.
#' @param H_strata a vector of population size by strata
#' @param J_strata a matrix whose columns sum to one giving the distribution of
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
