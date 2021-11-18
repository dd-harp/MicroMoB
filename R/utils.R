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
