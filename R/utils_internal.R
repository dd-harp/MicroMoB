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


#' @noRd
time_varying_parameter <- function(param, p, tmax) {
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
    stopifnot(length(param) == p)
    if (p > 1) {
      out <- replicate(n = tmax, expr = param)
    } else {
      out <- matrix(data = param, nrow = 1, ncol = tmax)
    }
  }
  stopifnot(nrow(out) == p)
  stopifnot(ncol(out) == tmax)
  return(out)
}
