#' @title Setup humans with SIS pathogen model
#' @description A simple SIS (Susceptible-Infected-Susceptible) model
#' @param stochastic should the model update deterministically or stochastically?
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @param theta a time spent matrix
#' @param wf biting weights
#' @param H vector of strata population sizes
#' @param X number of infectious persons in each strata
#' @param b transmission efficiency (mosquito to human)
#' @param c transmission efficiency (human to mosquito)
#' @param r recovery rate (inverse of infectious duration)
#' @return no return value
#' @export
setup_humans_SIS <- function(model, stochastic, theta, wf = NULL, H, X, b = 0.55, c = 0.15, r = 1/200) {
  stopifnot(inherits(model, "MicroMoB"))
  stopifnot(inherits(theta, "matrix"))

  stopifnot(nrow(theta) >= ncol(theta))
  stopifnot(approx_equal(rowSums(theta), 1))

  n <- nrow(theta)
  p <- ncol(theta)
  stopifnot(p == model$global$p)

  stopifnot(length(H) == n)
  stopifnot(length(X) == n)
  stopifnot(X <= H)

  stopifnot(is.finite(c(b, c, r)))
  stopifnot(c(b, c, r) >= 0)
  stopifnot(c(b, c) <= 1)

  model$global$n <- n

  if (is.null(wf)) {
    wf <- rep(1, n)
  }

  stopifnot(length(wf) == n)
  stopifnot(is.finite(wf))
  stopifnot(wf >= 0)

  human_class <- c("SIS")
  if (stochastic) {
    human_class <- c(human_class, "SIS_stochastic")
  } else {
    human_class <- c(human_class, "SIS_deterministic")
  }

  model$human <- structure(list(), class = human_class)
  model$human$theta <- theta
  model$human$wf <- wf
  model$human$H <- H
  model$human$X <- X

  model$human$EIR <- rep(0, n)

  model$human$b <- b
  model$human$c <- c
  model$human$r <- r
}


#' @title Get parameters for SIS human model
#' @description The JSON config file should have 8 entries:
#'  * stochastic: a boolean value
#'  * theta: matrix (row major)
#'  * wf: vector
#'  * H: vector
#'  * X: vector
#'  * b: scalar
#'  * c: scalar
#'  * r: scalar
#'
#' For interpretation of the entries, please read [MicroMoB::setup_humans_SIS].
#' @param path a file path to a JSON file
#' @return a named [list]
#' @importFrom jsonlite read_json
#' @examples
#' # to see an example of proper JSON input, run the following
#' library(jsonlite)
#' n <- 6 # number of human population strata
#' p <- 5 # number of patches
#' theta <- matrix(rexp(n*p), nrow = n, ncol = p)
#' theta <- theta / rowSums(theta)
#' H <- rep(10, n)
#' X <- rep(3, n)
#' par <- list(
#'  "stochastic" = FALSE,
#'  "theta" = theta,
#'  "wf" = rep(1, n),
#'  "H" = H,
#'  "X" = X,
#'  "b" = 0.55,
#'  "c" = 0.15,
#'  "r" = 1/200
#' )
#' toJSON(par, pretty = TRUE)
#' @export
get_config_humans_SIS <- function(path) {
  pars <- read_json(path = file.path(path), simplifyVector = TRUE)

  stopifnot(length(pars) == 8L)
  stopifnot(is.logical(pars$stochastic))
  stopifnot(length(pars$stochastic) == 1L)

  stopifnot(is.numeric(pars$theta))
  stopifnot(is.matrix(pars$theta))

  stopifnot(is.numeric(pars$wf))
  stopifnot(is.vector(pars$wf))

  stopifnot(is.numeric(pars$H))
  stopifnot(is.vector(pars$H))

  stopifnot(is.numeric(pars$X))
  stopifnot(is.vector(pars$X))

  stopifnot(is.numeric(pars$b))
  stopifnot(is.numeric(pars$c))
  stopifnot(is.numeric(pars$r))

  return(pars)
}


# step (update)

#' @title Update SIS human model
#' @inheritParams step_humans
#' @return no return value
#' @export
step_humans.SIS <- function(model) {
  NextMethod()
}

#' @title Update SIS human model (deterministic)
#' @inheritParams step_humans
#' @return no return value
#' @importFrom stats pexp
#' @export
step_humans.SIS_deterministic <- function(model) {

  h <- model$human$EIR * model$human$b

  new_infections <- pexp(q = h) * (model$human$H - model$human$X)
  old_infections <- (1 - pexp(q = model$human$r)) * model$human$X

  model$human$X <- new_infections + old_infections

}

#' @title Update SIS human model (stochastic)
#' @inheritParams step_humans
#' @return no return value
#' @importFrom stats pexp rbinom
#' @export
step_humans.SIS_stochastic <- function(model) {

  h <- model$human$EIR * model$human$b
  n <- model$global$n

  new_infections <- rbinom(n = n, prob = pexp(q = h), size = model$human$H - model$human$X)
  old_infections <- rbinom(n = n, prob = 1 - pexp(q = model$human$r), size = model$human$X)

  model$human$X <- new_infections + old_infections

}


#' @title Compute human biting weights for SIS model (\eqn{w_{f}})
#' @inheritParams compute_wf
#' @return a vector of length `n` giving the biting weights of human hosts in each stratum
#' @export
compute_wf.SIS <- function(model) {
  model$human$wf
}


#' @title Compute net infectiousness for SIS model (\eqn{x})
#' @inheritParams compute_x
#' @return a vector of length `n` giving the net infectiousness of human hosts in each stratum
#' @export
compute_x.SIS <- function(model) {
  x <- (model$human$X / model$human$H) * model$human$c
  return(as.vector(x))
}


#' @title Compute human population strata sizes for SIS model (\eqn{H})
#' @inheritParams compute_H
#' @return a vector of length `n` giving the size of each human population stratum
#' @export
compute_H.SIS <- function(model) {
  model$human$H
}


#' @title Compute time at risk matrix for SIS model (\eqn{\Psi})
#' @inheritParams compute_Psi
#' @return a matrix with `n` rows and `p` columns, the time at risk matrix
#' @export
compute_Psi.SIS <- function(model) {
  model$human$theta
}

