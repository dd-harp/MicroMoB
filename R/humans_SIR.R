#' @title Setup humans with SIR infection model
#' @description A simple SIR (Susceptible-Infected-Recovered) model
#' @param stochastic should the model update deterministically or stochastically?
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @param theta a time spent matrix
#' @param wf biting weights
#' @param H vector of strata population sizes
#' @param SIR a matrix giving S, I, R counts (columns) for each strata (rows)
#' @param b transmission efficiency (mosquito to human)
#' @param c transmission efficiency (human to mosquito)
#' @param gamma rate of recovery
#' @return no return value
#' @export
setup_humans_SIR <- function(model, stochastic, theta, wf = NULL, H, SIR, b = 0.55, c = 0.15, gamma = 1/5) {
  stopifnot(inherits(model, "MicroMoB"))
  stopifnot(inherits(theta, "matrix"))
  stopifnot(inherits(SIR, "matrix"))

  stopifnot(nrow(theta) >= ncol(theta))
  stopifnot(approx_equal(rowSums(theta), 1))

  n <- nrow(theta)
  p <- ncol(theta)
  stopifnot(p == model$global$p)

  stopifnot(length(H) == n)

  stopifnot(nrow(SIR) == n)
  stopifnot(ncol(SIR) == 3L)
  stopifnot(rowSums(SIR) == H)

  stopifnot(is.finite(c(b, c, gamma)))
  stopifnot(c(b, c, gamma) >= 0)
  stopifnot(c(b, c) <= 1)

  model$global$n <- n

  if (is.null(wf)) {
    wf <- rep(1, n)
  }

  stopifnot(length(wf) == n)
  stopifnot(is.finite(wf))
  stopifnot(wf >= 0)

  if (is.null(colnames(SIR))) {
    colnames(SIR) <- c("S", "I", "R")
  } else {
    stopifnot(colnames(SIR) == c("S", "I", "R"))
  }

  human_class <- c("SIR")
  if (stochastic) {
    human_class <- c(human_class, "SIR_stochastic")
  } else {
    human_class <- c(human_class, "SIR_deterministic")
  }

  model$human <- structure(list(), class = human_class)
  model$human$theta <- theta
  model$human$wf <- wf
  model$human$H <- H
  model$human$SIR <- SIR

  model$human$EIR <- rep(0, n)

  model$human$b <- b
  model$human$c <- c
  model$human$gamma <- gamma
}


#' @title Get parameters for SIR human model
#' @description The JSON config file should have 8 entries:
#'  * stochastic: a boolean value
#'  * theta: matrix (row major)
#'  * wf: vector
#'  * H: vector
#'  * SIR: matrix (row major)
#'  * b: scalar
#'  * c: scalar
#'  * gamma: scalar
#'
#' For interpretation of the entries, please read [MicroMoB::setup_humans_SIR].
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
#' SIR <- matrix(0, nrow = n, ncol = 3)
#' SIR[, 1] <- H
#' par <- list(
#'  "stochastic" = FALSE,
#'  "theta" = theta,
#'  "wf" = rep(1, n),
#'  "H" = H,
#'  "SIR" = SIR,
#'  "b" = 0.55,
#'  "c" = 0.15,
#'  "gamma" = 1/7
#' )
#' toJSON(par, pretty = TRUE)
#' @export
get_config_humans_SIR <- function(path) {
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

  stopifnot(is.numeric(pars$SIR))
  stopifnot(is.matrix(pars$SIR))

  stopifnot(is.numeric(pars$b))
  stopifnot(is.numeric(pars$c))
  stopifnot(is.numeric(pars$gamma))

  return(pars)
}



# step (update)

#' @title Update SIR human model
#' @inheritParams step_humans
#' @return no return value
#' @export
step_humans.SIR <- function(model) {
  NextMethod()
}

#' @title Update SIR human model (deterministic)
#' @inheritParams step_humans
#' @return no return value
#' @importFrom stats pexp
#' @export
step_humans.SIR_deterministic <- function(model) {

  h <- model$human$EIR * model$human$b

  # compute differences: S
  S_leave <- model$human$SIR[, "S"] * pexp(q =  h)

  # compute differences: I
  I_leave <- model$human$SIR[, "I"] * pexp(q =  model$human$gamma)

  # update
  model$human$SIR[, "S"] <- model$human$SIR[, "S"] - S_leave
  model$human$SIR[, "I"] <- model$human$SIR[, "I"] + S_leave - I_leave
  model$human$SIR[, "R"] <- model$human$SIR[, "R"] + I_leave

  model$human$SIR <- pmax(model$human$SIR, 0)

}

#' @title Update SIR human model (stochastic)
#' @inheritParams step_humans
#' @return no return value
#' @importFrom stats pexp rbinom
#' @export
step_humans.SIR_stochastic <- function(model) {

  h <- model$human$EIR * model$human$b

  n <- model$global$n

  # compute differences: S
  S_leave <- rbinom(n = n, size = model$human$SIR[, "S"], prob = pexp(q =  h))

  # compute differences: I
  I_leave <- rbinom(n = n, size = model$human$SIR[, "I"], prob = pexp(q =  model$human$gamma))

  # update
  model$human$SIR[, "S"] <- model$human$SIR[, "S"] - S_leave
  model$human$SIR[, "I"] <- model$human$SIR[, "I"] + S_leave - I_leave
  model$human$SIR[, "R"] <- model$human$SIR[, "R"] + I_leave

}


#' @title Compute human biting weights for SIR model (\eqn{w_{f}})
#' @inheritParams compute_wf
#' @return a vector of length `n` giving the biting weights of human hosts in each stratum
#' @export
compute_wf.SIR <- function(model) {
  model$human$wf
}

#' @title Compute net infectiousness for SIR model (\eqn{x})
#' @inheritParams compute_x
#' @return a vector of length `n` giving the net infectiousness of human hosts in each stratum
#' @export
compute_x.SIR <- function(model) {
  X <- model$human$SIR[, "I"] / model$human$H
  X[is.nan(X)] <- 0
  return(as.vector(X * model$human$c))
}

#' @title Compute human population strata sizes for SIR model (\eqn{H})
#' @inheritParams compute_H
#' @return a vector of length `n` giving the size of each human population stratum
#' @export
compute_H.SIR <- function(model) {
  model$human$H
}


#' @title Compute time at risk matrix for SIR model (\eqn{\Psi})
#' @inheritParams compute_Psi
#' @return a matrix with `n` rows and `p` columns, the time at risk matrix
#' @export
compute_Psi.SIR <- function(model) {
  model$human$theta
}
