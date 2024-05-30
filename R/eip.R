# generic methods to compute the extrinsic incubation period (EIP)

#' @title Compute the EIP
#' @description This method dispatches on the type of `EIPmod`.
#' @param t current simulation time
#' @param MYZpar a [list]
#' @return [numeric]
#' @export
EIP <- function(t, MYZpar) {
  UseMethod("EIP", MYZpar$EIPmod)
}

#' @title Set up the fixed model for control forcing (do nothing)
#' @param EIPname the class name of the function
#' @param MYZpar the MYZ parameters
#' @param MYZopts is a [list] that overwrites default options
#' @return [list] MYZpar with the EIPmod attached
#' @export
setup_EIP <- function(EIPname, MYZpar, MYZopts = list()) {
  class(EIPname) <- EIPname
  UseMethod("setup_EIP", EIPname)
}

