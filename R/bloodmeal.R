#' @title Compute bloodmeals taken by mosquitoes on hosts
#' @description This should be run prior to any `step` functions to update
#' components over a time step. It computes various quantities related to
#' disease transmission between species using the generic interfaces (methods)
#' provided by each component. It updates the `EIR` vector for the human component, and `kappa`, the net infectiousness
#' of hosts for the mosquito component.
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @export
compute_bloodmeal <- function(model) {
  stopifnot(inherits(model, "MicroMoB"))

  # human quantities
  W <- compute_W(model)
  H <- compute_H(model)
  x <- compute_x(model)
  wf <- compute_wf(model)
  Psi <- compute_Psi(model)

  # biting distribution matrix (n x p)
  beta <- diag(wf) %*% Psi %*% diag(1/W)

  stopifnot(nrow(beta) == model$global$n)
  stopifnot(ncol(beta) == model$global$p)

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
  q <- (W + Wd) / B

  # density of infective mosquitoes
  Z <- compute_Z(model)

  # calculate EIR and kappa (mosy->human, human->mosy
  model$human$EIR <- beta %*% (f*q*v*Z)
  model$mosquito$kappa <- (v * t(beta) %*% (x*H)) + ((1 - v) * xd)

  stopifnot(length(model$human$EIR) == model$global$n)
  stopifnot(length(model$mosquito$kappa) == model$global$p)

}
