# utilities that are not documented

#' @title Does a numeric object consist of only zeros and ones?
#' @param x a [numeric] object
#' @return a logical value
#' @export
is_binary <- function(x) {
  uniq <- unique(as.vector(x))
  if (all(length(uniq) == 2L)) {
    return(all(uniq %in% c(0, 1)))
  } else {
    return(FALSE)
  }
}


#' @title Check if two numeric values are approximately equal
#' @param a a [numeric] object
#' @param b a [numeric] object
#' @param tol the numeric tolerance
#' @return a logical value
#' @export
approx_equal <- function(a, b, tol = sqrt(.Machine$double.eps)) {
  abs(a - b) < tol
}


#' @title Division of integers
#' @param a the dividend
#' @param b the divisor
#' @return a list with two elements, `quo` (quotient) and `rem` (remainder)
#' @export
divmod <- function(a,b){
  a <- as.integer(a)
  b <- as.integer(b)
  c(
    quo = a %/% b,
    rem = a %% b
  )
}


#' @title Distribute items into bins as evenly as possible
#' @param n number of bins
#' @param p number of items
#' @return a numeric vector of bin sizes
#' @export
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
