# generic methods for adult component

#' @title Set bloodfeeding and mortality rates to baseline
#' @description This method dispatches on the type of `pars$MYZpar`. It should
#' set the values of the bionomic parameters to baseline values.
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @param s the species index
#' @return a [list]
#' @export
MBionomics <- function(t, y, pars, s) {
  UseMethod("MBionomics", pars$MYZpar[[s]])
}

#' @title Time spent host seeking/feeding and resting/ovipositing
#' @description This method dispatches on the type of `pars$MYZpar`.
#' @param t current simulation time
#' @param pars a [list]
#' @return either a [numeric] vector if the model supports this feature, or [NULL]
#' @export
F_tau <- function(t, pars) {
  UseMethod("F_tau", pars$MYZpar)
}

#' @title Blood feeding rate of the infective mosquito population
#' @description This method dispatches on the type of `pars$MYZpar`.
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @param s the species index
#' @return a [numeric] vector of length `nPatches`
#' @export
F_fqZ <- function(t, y, pars, s) {
  UseMethod("F_fqZ", pars$MYZpar[[s]])
}

#' @title Blood feeding rate of the mosquito population
#' @description This method dispatches on the type of `pars$MYZpar`.
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @param s the species index
#' @return a [numeric] vector of length `nPatches`
#' @export
F_fqM <- function(t, y, pars, s) {
  UseMethod("F_fqM", pars$MYZpar[[s]])
}

#' @title Number of eggs laid by adult mosquitoes
#' @description This method dispatches on the type of `pars$MYZpar`.
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @param s the species index
#' @return a [numeric] vector of length `nPatches`
#' @export
F_eggs <- function(t, y, pars, s) {
  UseMethod("F_eggs", pars$MYZpar[[s]])
}


#' @title Derivatives for adult mosquitoes
#' @description This method dispatches on the type of `pars$MYZpar`.
#' @param t current simulation time
#' @param y state vector
#' @param pars a [list]
#' @param s the species index
#' @return the derivatives a [vector]
#' @export
dMYZdt <- function(t, y, pars, s) {
  UseMethod("dMYZdt", pars$MYZpar[[s]])
}

#' @title A function to set up adult mosquito models
#' @description This method dispatches on `MYZname`.
#' @param MYZname the name of the model
#' @param pars a [list]
#' @param s the species index
#' @param MYZopts a [list]
#' @param EIPname is the name of the EIPmod
#' @param calK is a [matrix]
#' @return [list]
#' @export
setup_MYZpar = function(MYZname, pars, s, MYZopts=list(), EIPname="fixed", calK=diag(1)){
  class(MYZname) <- MYZname
  UseMethod("setup_MYZpar", MYZname)
}

#' @title A function to set up adult mosquito models
#' @description This method dispatches on `MYZname`.
#' @param pars a [list]
#' @param s the species index
#' @param MYZopts a [list]
#' @return [list]
#' @export
setup_MYZinits = function(pars, s, MYZopts=list()){
  UseMethod("setup_MYZinits", pars$MYZpar[[s]])
}

#' @title Return the variables as a list
#' @description This method dispatches on the type of `pars$MYZpar[[s]]`.
#' @param y the variables
#' @param pars a [list]
#' @param s the vector species index
#' @return a [list]
#' @export
list_MYZvars <- function(y, pars, s) {
  UseMethod("list_MYZvars", pars$MYZpar[[s]])
}

#' @title Add indices for adult mosquitoes to parameter list
#' @description This method dispatches on the type of `pars$MYZpar`.
#' @param pars a [list]
#' @param s the species index
#' @return [list]
#' @export
make_indices_MYZ <- function(pars, s) {
  UseMethod("make_indices_MYZ", pars$MYZpar[[s]])
}

#' @title Parse the output of deSolve and return the variables by name in a list
#' @description This method dispatches on the type of `pars$MYZpar`.
#' It computes the variables by name and returns a named list.
#' @param dts_out a [matrix] of outputs from deSolve
#' @param pars a [list] that defines a model
#' @param s the species index
#' @return [list]
#' @export
parse_dts_out_MYZ <- function(dts_out, pars, s) {
  UseMethod("parse_dts_out_MYZ", pars$MYZpar[[s]])
}

#' @title Return initial values as a vector
#' @description This method dispatches on the type of `pars$MYZpar`.
#' @param pars a [list]
#' @param s the species index
#' @return [numeric]
#' @export
get_inits_MYZ <- function(pars, s) {
  UseMethod("get_inits_MYZ", pars$MYZpar[[s]])
}

#' @title Set the initial values as a vector
#' @description This method dispatches on the type of `pars$MYZpar`.
#' @param pars a [list]
#' @param y0 a vector of variable values from a simulation
#' @param s the species index
#' @return a [list]
#' @export
update_inits_MYZ <- function(pars, y0, s) {
  UseMethod("update_inits_MYZ", pars$MYZpar[[s]])
}

#' @title Make the mosquito demography matrix
#' @param t the current time
#' @param MYZpar a list defining a mosquito model
#' @return a [matrix] of dimensions `nPatches` by `nPatches`
#' @export
make_Omega <- function(t, MYZpar) {with(MYZpar,
  diag(p*(1-sigma), nPatches) + diag(p*sigma, nPatches) %*% calK
)}

