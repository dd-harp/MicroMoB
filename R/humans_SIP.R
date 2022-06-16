#' @title Setup humans with SIP pathogen model
#' @description A simple SIP (Susceptible-Infected-Protected) model
#' @param stochastic should the model update deterministically or stochastically?
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @param theta a time spent matrix
#' @param wf biting weights
#' @param SIP matrix of strata (rows) by health states (SIP)
#' @param b transmission efficiency (mosquito to human)
#' @param c transmission efficiency (human to mosquito)
#' @param r recovery rate (inverse of infectious duration)
#' @param rho probability of treatment upon infection
#' @param eta rate at which prophylaxis decays
#' @return no return value
#' @export
setup_humans_SIP <- function(model, stochastic, theta, wf = NULL, SIP, b = 0.55, c = 0.15, r = 1/200, rho = 0.07, eta = 1/32) {
  stopifnot(inherits(model, "MicroMoB"))
  stopifnot(inherits(theta, "matrix"))

  stopifnot(nrow(theta) >= ncol(theta))
  stopifnot(approx_equal(rowSums(theta), 1))

  n <- nrow(theta)
  p <- ncol(theta)
  stopifnot(p == model$global$p)

  stopifnot(nrow(SIP) == n)
  stopifnot(ncol(SIP) == 3L)
  stopifnot(SIP >= 0)
  colnames(SIP) <- c("S", "I", "P")

  stopifnot(is.finite(c(b, c, r, rho, eta)))
  stopifnot(c(b, c, r, rho, eta) >= 0)
  stopifnot(c(b, c, rho, eta) <= 1)

  model$global$n <- n

  if (is.null(wf)) {
    wf <- rep(1, n)
  }

  stopifnot(length(wf) == n)
  stopifnot(is.finite(wf))
  stopifnot(wf >= 0)

  human_class <- c("SIP")
  if (stochastic) {
    human_class <- c(human_class, "SIP_stochastic")
  } else {
    human_class <- c(human_class, "SIP_deterministic")
  }

  model$human <- structure(list(), class = human_class)
  model$human$theta <- theta
  model$human$wf <- wf
  model$human$SIP <- SIP

  model$human$EIR <- rep(0, n)
  model$human$incidence <- rep(0, n)

  model$human$b <- b
  model$human$c <- c
  model$human$r <- r
  model$human$rho <- rho
  model$human$eta <- eta
}


# step (update)

#' @title Update SIP human model
#' @inheritParams step_humans
#' @return no return value
#' @export
step_humans.SIP <- function(model) {
  NextMethod()
}

#' @title Update SIP human model (deterministic)
#' @inheritParams step_humans
#' @return no return value
#' @importFrom stats pexp
#' @export
step_humans.SIP_deterministic <- function(model) {

  h <- model$human$EIR * model$human$b

  new_infections <- pexp(q = h) * model$human$SIP[, "S"]
  new_treatment <- new_infections * model$human$rho
  new_disease <- new_infections - new_treatment

  new_recoveries <- pexp(q = model$human$r) * model$human$SIP[, "I"]

  new_wane <- pexp(q = model$human$eta) * model$human$SIP[, "P"]

  model$human$SIP[, "S"] <- model$human$SIP[, "S"] - new_disease - new_treatment + new_recoveries + new_wane
  model$human$SIP[, "I"] <- model$human$SIP[, "I"] + new_disease - new_recoveries
  model$human$SIP[, "P"] <- model$human$SIP[, "P"] + new_treatment - new_wane

  model$human$incidence <- new_infections

}

#' @title Update SIP human model (stochastic)
#' @inheritParams step_humans
#' @return no return value
#' @importFrom stats pexp rbinom
#' @export
step_humans.SIP_stochastic <- function(model) {

  h <- model$human$EIR * model$human$b
  n <- model$global$n

  new_infections <- rbinom(n = n, prob = pexp(q = h), size = model$human$SIP[, "S"])
  new_treatment <- rbinom(n = n, size = new_infections, prob = pexp(q = model$human$rho))
  new_disease <- new_infections - new_treatment

  new_recoveries <- rbinom(n = n, size = model$human$SIP[, "I"], prob = pexp(q = model$human$r))

  new_wane <- rbinom(n = n, size = model$human$SIP[, "P"], prob = pexp(q = model$human$eta))

  model$human$SIP[, "S"] <- model$human$SIP[, "S"] - new_disease - new_treatment + new_recoveries + new_wane
  model$human$SIP[, "I"] <- model$human$SIP[, "I"] + new_disease - new_recoveries
  model$human$SIP[, "P"] <- model$human$SIP[, "P"] + new_treatment - new_wane

  model$human$incidence <- new_infections

}


#' @title Compute human biting weights for SIP model (\eqn{w_{f}})
#' @inheritParams compute_wf
#' @return a vector of length `n` giving the biting weights of human hosts in each stratum
#' @export
compute_wf.SIP <- function(model) {
  model$human$wf
}


#' @title Compute net infectiousness for SIP model (\eqn{x})
#' @inheritParams compute_x
#' @return a vector of length `n` giving the net infectiousness of human hosts in each stratum
#' @export
compute_x.SIP <- function(model) {
  x <- (model$human$SIP[2L, ] / rowSums(model$human$SIP)) * model$human$c
  return(as.vector(x))
}


#' @title Compute human population strata sizes for SIP model (\eqn{H})
#' @inheritParams compute_H
#' @return a vector of length `n` giving the size of each human population stratum
#' @export
compute_H.SIP <- function(model) {
  rowSums(model$human$SIP)
}


#' @title Compute time at risk matrix for SIP model (\eqn{\Psi})
#' @inheritParams compute_Psi
#' @return a matrix with `n` rows and `p` columns, the time at risk matrix
#' @export
compute_Psi.SIP <- function(model) {
  model$human$theta
}

