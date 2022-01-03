# utilities that are not documented or exported

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


#' @title Input parameters that may vary by time and patch
#' @param param if given a matrix, it must have nrows equal to `p` and ncols equal to
#' either `tmax` or `365`; if given a vector it must be of length `p`, `tmax`, or `365`.
#' @param p number of patches
#' @param tmax number of time steps
time_patch_varying_parameter <- function(param, p, tmax) {
  if (inherits(param, "matrix")) {
    stopifnot(nrow(param) == p)
    if (ncol(param) == 365L) {
      ix <- (1:tmax) %% 365L
      ix[which(ix == 0L)] <- 365L
      out <- param[, ix, drop = FALSE]
    } else if (ncol(param) == tmax) {
      out <- param
    } else {
      stop("incorrect dimensions of parameter")
    }
  } else {
    # vector input
    if (length(param) == p) {
      if (p > 1) {
        out <- replicate(n = tmax, expr = param)
      } else {
        out <- matrix(data = param, nrow = 1, ncol = tmax)
      }
    } else if (length(param) == tmax) {
      out <- do.call(rbind, replicate(n = p, expr = param, simplify = FALSE))
    } else if (length(param) == 365L) {
      ix <- (1:tmax) %% 365L
      ix[which(ix == 0L)] <- 365L
      out <- do.call(rbind, replicate(n = p, expr = param[ix], simplify = FALSE))
    } else {
      stop("incorrect length of parameter")
    }
  }
  stopifnot(nrow(out) == p)
  stopifnot(ncol(out) == tmax)
  return(out)
}


#' @title Input parameters that may vary by time
#' @param param a vector of length `1`, `tmax`, or `365`.
#' @param tmax number of time steps
time_varying_parameter <- function(param, tmax) {
  if (length(param) == 1L) {
    out <- rep(param, tmax)
  } else if(length(param) == 365L) {
    ix <- (1:tmax) %% 365L
    ix[which(ix == 0L)] <- 365L
    out <- param[ix]
  } else if(length(param) == tmax) {
    out <- param
  } else {
    stop("incorrect length of parameter")
  }
  stopifnot(length(out) == tmax)
  return(out)
}
