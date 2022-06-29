#' @title Compute bloodmeals taken by mosquitoes on hosts
#' @description This should be run prior to any `step` functions to update
#' components over a time step. It computes various quantities related to
#' disease transmission between species using the generic interfaces (methods)
#' provided by each component. It updates the `EIR` vector for the human component, and `kappa`, the net infectiousness
#' of hosts for the mosquito component.
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @return no return value
#' @export
compute_bloodmeal <- function(model) {
  stopifnot(inherits(model, "MicroMoB"))

  n <- model$global$n
  p <- model$global$p

  # human quantities
  H <- compute_H(model)
  x <- compute_x(model)
  wf <- compute_wf(model)
  Psi <- compute_Psi(model)
  W <- as.vector(t(Psi) %*% (wf * H))

  # biting distribution matrix (n x p)
  beta <- diag(wf, nrow = n, ncol = n) %*% Psi %*% diag(1/W, nrow = p, ncol = p)

  stopifnot(nrow(beta) == n)
  stopifnot(ncol(beta) == p)

  # host availability
  Wd <- compute_Wd(model)
  O <- compute_O(model)
  B <- W + Wd + O

  # visitor net infectiousness
  xd <- compute_xd(model)

  # blood feeding rates
  f <- compute_f(model, B = B)

  # fraction feeding on resident humans
  v <- W / (W + Wd)

  # human blood feeding fraction
  q <- compute_q(model, W = W, Wd = Wd, B = B)

  # density of infective mosquitoes
  Z <- compute_Z(model)

  # calculate EIR and kappa (mosy->human, human->mosy
  model$human$EIR <- beta %*% (f*q*v*Z)
  model$human$EIR <- as.vector(model$human$EIR)

  model$mosquito$kappa <- (v * (t(beta) %*% (x*H))) + ((1 - v) * xd)
  model$mosquito$kappa <- as.vector(model$mosquito$kappa)

  # mosquito feeding habit: f, q
  model$mosquito$f <- f
  model$mosquito$q <- q

  # check dimensions
  stopifnot(length(model$human$EIR) == n)
  stopifnot(length(model$mosquito$kappa) == p)
  stopifnot(length(f) == p)
  stopifnot(length(q) == p)

}


#' @title Compute bloodmeals taken by mosquitoes on hosts in simple models
#' @description The difference between this and [compute_bloodmeal] is that
#' this function does not include any computations of alterative blood hosts
#' or visitors and is suitable for models which only include mosquitoes and resident
#' human populations.
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @return no return value
#' @export
compute_bloodmeal_simple <- function(model) {
  stopifnot(inherits(model, "MicroMoB"))

  n <- model$global$n
  p <- model$global$p

  # human quantities
  H <- compute_H(model)
  x <- compute_x(model)
  wf <- compute_wf(model)
  Psi <- compute_Psi(model)
  W <- as.vector(t(Psi) %*% (wf * H))

  # biting distribution matrix (n x p)
  beta <- diag(wf, nrow = n, ncol = n) %*% Psi %*% diag(1/W, nrow = p, ncol = p)

  stopifnot(nrow(beta) == n)
  stopifnot(ncol(beta) == p)
  empty_vec <- rep(0, p)

  # blood feeding rates
  f <- compute_f(model, B = empty_vec)

  # human blood feeding fraction
  q <- compute_q(model, W = W, Wd = empty_vec, B = empty_vec)

  # density of infective mosquitoes
  Z <- compute_Z(model)

  # calculate EIR and kappa (mosy->human, human->mosy
  model$human$EIR <- beta %*% (f*q*Z)
  model$human$EIR <- as.vector(model$human$EIR)

  model$mosquito$kappa <- t(beta) %*% (x*H)
  model$mosquito$kappa <- as.vector(model$mosquito$kappa)

  # mosquito feeding habit: f, q
  model$mosquito$f <- f
  model$mosquito$q <- q

  # check dimensions
  stopifnot(length(model$human$EIR) == n)
  stopifnot(length(model$mosquito$kappa) == p)
  stopifnot(length(f) == p)
  stopifnot(length(q) == p)

}
