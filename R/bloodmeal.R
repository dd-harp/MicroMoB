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


}
