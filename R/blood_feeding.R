


#' @title Compute biting distribution matrix (\eqn{\beta})
#' @param human an object from [MicroMoB::setup_human]
#' @param xi probabilities to initiate a feeding search during each fraction
#' of the day (must sum to 1)
#' @param t time
#' @seealso [MicroMoB::compute_beta.day] [MicroMoB::compute_beta.dt]
#' @export
compute_beta <- function(human, xi, t) {
  UseMethod("compute_beta", human$timespent)
}

#' @title Compute daily biting distribution matrix (\eqn{\beta})
#' @inheritParams compute_beta
#' @return a matrix of dimension \eqn{n \times p}
#' @export
compute_beta.day <- function(human, xi = 1, t) {
  Psi <- t(compute_Psi(human = human, xi = xi, t = t))
  wf <- compute_wf(biteweight = human$biteweight, t = t)
  W <- compute_W(human = human, Psi_t = Psi_t, t = t)
  return(diag(wf) %*% Psi %*% diag(1/W))
}

#' @title Compute fractional biting distribution array (\eqn{\beta})
#' @inheritParams compute_beta
#' @return an array of dimension \eqn{n \times p \times d}
#' @export
compute_beta.dt <- function(human, xi, t) {

  beta <- array(data = 0, dim = c(human$n, human$p, human$timespent$d))

  Psi_t <- compute_Psi(human = human, xi = xi, t = t)
  wf <- compute_wf(biteweight = human$biteweight, t = t)
  wf <- diag(wf)
  W <- compute_W(human = human, Psi_t = Psi_t, t = t)

  for (k in 1:human$timespent$d) {
    beta[, , k] <- wf %*% t(Psi_t[, , k]) %*% diag(W[, k])
  }

  return(beta)
}

#' @export
compute_beta.default <- function(human, xi, t) {
  stop("compute_beta has no method for dispatch type ", class(human$timespent))
}
