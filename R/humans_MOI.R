# a MoI (multiplicity of infection) model for the human component

#' @title Setup humans with MOI (multiplicity of infection) pathogen model
#' @description This is a queueing model (M/M/inf) of superinfection in humans.
#' @note The [MicroMoB::step_humans] method for the MOI model will grow the `MOI`
#' matrix (add rows) if an individual's MOI exceeds the size of the matrix; therefore
#' it's a good idea to pad the input matrix with extra empty rows to avoid
#' reallocating memory during the simulation as much as possible.
#' @param stochastic should the model update deterministically or stochastically?
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @param theta a time spent matrix
#' @param wf biting weights
#' @param H vector of strata population sizes
#' @param MOI a matrix giving the distribution of persons across strata (columns) and
#' multiplicity of infection (rows).
#' @param b transmission efficiency (mosquito to human)
#' @param c transmission efficiency (human to mosquito)
#' @param r recovery rate (inverse of infectious duration)
#' @param sigma control non-independence of pathogen clearance; `sigma > 1` indicates competition
#' (clearance is faster than independent) and `sigma < 1` indicates facilitation (clearance is slower than independent).
#' @return no return value
#' @export
setup_humans_MOI <- function(model, stochastic, theta, wf = NULL, H, MOI, b = 0.55, c = 0.15, r = 1/200, sigma = 1) {
  stopifnot(inherits(model, "MicroMoB"))
  stopifnot(inherits(theta, "matrix"))

  stopifnot(nrow(theta) >= ncol(theta))
  stopifnot(approx_equal(rowSums(theta), 1))

  n <- nrow(theta)
  p <- ncol(theta)
  stopifnot(p == model$global$p)

  stopifnot(length(H) == n)

  stopifnot(colSums(MOI) == H)
  stopifnot(ncol(MOI) == n)
  stopifnot(nrow(MOI) > 2)

  stopifnot(is.finite(c(b, c, r, sigma)))
  stopifnot(c(b, c, r, sigma) >= 0)
  stopifnot(c(b, c) <= 1)

  model$global$n <- n

  if (is.null(wf)) {
    wf <- rep(1, n)
  }

  stopifnot(length(wf) == n)
  stopifnot(is.finite(wf))
  stopifnot(wf >= 0)

  human_class <- c("MOI")
  if (stochastic) {
    human_class <- c(human_class, "MOI_stochastic")
  } else {
    human_class <- c(human_class, "MOI_deterministic")
  }

  model$human <- structure(list(), class = human_class)
  model$human$theta <- theta
  model$human$wf <- wf
  model$human$H <- H
  model$human$MOI <- MOI

  model$human$EIR <- rep(0, n)

  model$human$b <- b
  model$human$c <- c
  model$human$r <- r
  model$human$sigma <- sigma
}


# step (update)

#' @title Update MOI human model
#' @inheritParams step_humans
#' @return no return value
#' @export
step_humans.MOI <- function(model) {
  NextMethod()
}

#' @title Update MOI human model (deterministic)
#' @inheritParams step_humans
#' @return no return value
#' @importFrom stats pexp rbinom
#' @importFrom utils tail
#' @export
step_humans.MOI_deterministic <- function(model) {

  maxMOI <- nrow(model$human$MOI)

  h <- model$human$EIR * model$human$b # FOI by strata
  h <- pexp(q = h)

  rho <- model$human$r * (1:(maxMOI-1))^model$human$sigma # recovery by MOI
  rho <- pexp(q = rho)

  # new infections in each bin
  new_infections <- model$human$MOI %*% diag(h)

  # new recoveries in each bin
  recoveries <- diag(rho) %*% model$human$MOI[-1L, ]

  # check if we need to extend MOI
  tot_infections <- rowSums(new_infections)
  if (tail(tot_infections, 1L) > 0) {
    recoveries <- rbind(recoveries, 0)
    new_infections <- rbind(new_infections, 0)
    model$human$MOI <- rbind(model$human$MOI, 0)
    maxMOI <- nrow(model$human$MOI)
  }

  # apply state updates
  model$human$MOI <- model$human$MOI - new_infections
  model$human$MOI[-1L, ] <- model$human$MOI[-1L, ] + new_infections[-maxMOI, ]

  model$human$MOI[-1L, ] <- model$human$MOI[-1L, ] - recoveries
  model$human$MOI[-maxMOI, ] <- model$human$MOI[-maxMOI, ] + recoveries

}

#' @title Update MOI human model (stochastic)
#' @inheritParams step_humans
#' @return no return value
#' @importFrom stats pexp rbinom
#' @importFrom abind abind
#' @export
step_humans.MOI_stochastic <- function(model) {

  n <- model$global$n

  maxMOI <- nrow(model$human$MOI)

  h <- model$human$EIR * model$human$b # FOI by strata

  rho <- model$human$r * (1:(maxMOI-1))^model$human$sigma # recovery by MOI

  # who experiences some events in each strata
  events <- vapply(X = 1:n, FUN = function(i) {

    # P(infection or recovery)
    hazards <- rep(h[i], maxMOI)
    hazards[2:maxMOI] <- hazards[2:maxMOI] + rho
    probs <- pexp(q = hazards) # P(inf or rec) = 1 - exp(-inf + rec)

    # sample who experiences either event
    any <- rbinom(n = maxMOI, size = model$human$MOI[, i], prob = probs)

    # sample recovery (rho / rho + h)
    recovery <- rbinom(n = maxMOI, size = any, prob = c(0, rho) / hazards)
    infection <- any - recovery

    cbind(infection, recovery)

  }, FUN.VALUE = matrix(0, maxMOI, 2), USE.NAMES = FALSE)

  # check if we need to extend MOI
  if (any(events[maxMOI, 2, ] > 0)) {
    events <- abind(events, matrix(data = 0, nrow = 2, ncol = n), along = 1)
    model$human$MOI <- rbind(model$human$MOI, 0)
    maxMOI <- nrow(model$human$MOI)
  }

  # apply state updates
  for (i in 1:n) {
    # infections
    model$human$MOI[, i] <- model$human$MOI[, i] - events[, 1L, i]
    model$human$MOI[-1L, i] <- model$human$MOI[-1L, i] + events[-maxMOI, 1L, i]
    # recoveries
    model$human$MOI[-1L, i] <- model$human$MOI[-1L, i] - events[-1L, 2L, i]
    model$human$MOI[-maxMOI, i] <- model$human$MOI[-maxMOI, i] + events[-1L, 2L, i]
  }

}


# biting computation

#' @title Compute human biting weights for MOI model (\eqn{w_{f}})
#' @inheritParams compute_wf
#' @return a vector of length `n` giving the biting weights of human hosts in each stratum
#' @export
compute_wf.MOI <- function(model) {
  model$human$wf
}

#' @title Compute net infectiousness for MOI model (\eqn{x})
#' @description In the simple MOI (queueing) model here (M/M/inf), net infectiousness
#' is considered not to vary with increasing MOI. It is calculated as
#' \deqn{c \cdot (1 - \frac{X_{0}}{H})}
#' where \eqn{X_{0}} is the number of uninfected persons (multiplicity of infection of zero).
#' @inheritParams compute_x
#' @return a vector of length `n` giving the net infectiousness of human hosts in each stratum
#' @export
compute_x.MOI <- function(model) {
  X <- (model$human$H - model$human$MOI[1, ]) / model$human$H
  return(as.vector(X * model$human$c))
}

#' @title Compute human population strata sizes for MOI model (\eqn{H})
#' @inheritParams compute_H
#' @return a vector of length `n` giving the size of each human population stratum
#' @export
compute_H.MOI <- function(model) {
  model$human$H
}

#' @title Compute time at risk matrix for MOI model (\eqn{\Psi})
#' @inheritParams compute_Psi
#' @return a matrix with `n` rows and `p` columns, the time at risk matrix
#' @export
compute_Psi.MOI <- function(model) {
  model$human$theta
}

