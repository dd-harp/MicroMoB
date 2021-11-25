# beta: biting distribution matrix

#' @title Compute biting distribution matrix (\eqn{\beta})
#' @param human an object from [MicroMoB::setup_human]
#' @param Psi_t a time at risk matrix from [MicroMoB::compute_Psi]
#' @param t time
#' @seealso [MicroMoB::compute_beta.day] [MicroMoB::compute_beta.dt]
#' @export
compute_beta <- function(human, Psi_t, t) {
  UseMethod("compute_beta", human$timespent)
}

#' @title Compute daily biting distribution matrix (\eqn{\beta})
#' @inheritParams compute_beta
#' @return a matrix of dimension \eqn{n \times p}
#' @export
compute_beta.day <- function(human, Psi_t, t) {
  wf <- compute_wf(biteweight = human$biteweight, t = t)
  W <- compute_W(human = human, Psi_t = Psi_t, t = t)
  beta <- diag(wf) %*% t(Psi_t) %*% diag(1/W)
  return(beta)
}

#' @title Compute fractional biting distribution array (\eqn{\beta})
#' @inheritParams compute_beta
#' @return an array of dimension \eqn{n \times p \times d}
#' @export
compute_beta.dt <- function(human, Psi_t, t) {

  beta <- array(data = 0, dim = c(human$n, human$p, human$timespent$d))

  wf <- compute_wf(biteweight = human$biteweight, t = t)
  wf <- diag(wf)
  W <- compute_W(human = human, Psi_t = Psi_t, t = t)

  for (k in 1:human$timespent$d) {
    beta[, , k] <- wf %*% t(Psi_t[, , k]) %*% diag(1/W[, k])
  }

  return(beta)
}

#' @export
compute_beta.default <- function(human, Psi_t, t) {
  stop("compute_beta has no method for dispatch type ", class(human$timespent))
}


# parasite transmission: dispatch on the human model, which may be dt or day

#' @title Compute parasite transmission terms from bloodmeals (\eqn{\kappa}, \eqn{h})
#' @param model a model object (an [environment])
#' @param t time
#' @seealso [MicroMoB::compute_bloodmeal.dt] [MicroMoB::compute_bloodmeal.day]
#' @export
compute_bloodmeal <- function(model, t) {
  UseMethod("compute_bloodmeal", model$human$timespent)
}

#' @title Compute parasite transmission terms from bloodmeals with fractional human time spent (\eqn{\kappa}, \eqn{h})
#' @inheritParams compute_bloodmeal
#' @export
compute_bloodmeal.dt <- function(model, t) {

  xi <- model$mosquito$xi
  Psi_t <- compute_Psi(human = model$human, xi = xi, t = t)

  beta_array <- compute_beta(human = model$human, Psi_t = Psi_t, t = t)

  W <- compute_W(human = model$human, Psi_t = Psi_t, t = t)
  B <- compute_B(otherhosts = model$otherhosts, t = t)
  W_delta <- compute_W_delta(visitors = model$visitors, t = t)

  Tot <- W + W_delta + B
  upsilon <- W / W + W_delta

  f <- Ff_0(x = Tot, fx = model$bm$fx, sf = model$bm$sf)
  q <- W + W_delta / Tot


  # h=bβ·fqυZ

  # κ=υβT·x+ (1−υ)xδ

}

#' @title Compute parasite transmission terms from bloodmeals with daily human time spent (\eqn{\kappa}, \eqn{h})
#' @inheritParams compute_bloodmeal
#' @export
compute_bloodmeal.day <- function(model, t) {

}

#' @export
compute_bloodmeal.default <- function(model, t) {
  stop("compute_bloodmeal has no method for dispatch type ", class(model$human$timespent))
}
