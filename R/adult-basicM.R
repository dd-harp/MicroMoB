# specialized methods for the adult mosquito basicM model

#' @title Reset bloodfeeding and mortality rates to baseline
#' @description Implements [MBionomics] for the basicM model
#' @inheritParams MBionomics
#' @return the model as a [list]
#' @export
MBionomics.basicM <- function(t, y, pars, s) {
  with(pars$MYZpar[[s]],{
    pars$MYZpar[[s]]$f     <- F_f(t, pars$MYZpar[[s]])
    pars$MYZpar[[s]]$q     <- F_q(t, pars$MYZpar[[s]])
    pars$MYZpar[[s]]$p     <- F_p(t, pars$MYZpar[[s]])
    pars$MYZpar[[s]]$sigma <- F_sigma(t, pars$MYZpar[[s]])
    pars$MYZpar[[s]]$nu    <- F_nu(t, pars$MYZpar[[s]])
    pars$MYZpar[[s]]$Omega <- make_Omega(t, pars$MYZpar[[s]])
    return(pars)
  })}

#' @title The net blood feeding rate of the infective mosquito population in a patch
#' @description Implements [F_fqZ] for the basicM model.
#' @inheritParams F_fqZ
#' @return a [numeric] vector of length `nPatches`
#' @export
F_fqZ.basicM <- function(t, y, pars, s) {
  numeric(0)
}

#' @title The net blood feeding rate of the infective mosquito population in a patch
#' @description Implements [F_fqM] for the basicM model.
#' @inheritParams F_fqM
#' @return a [numeric] vector of length `nPatches`
#' @export
F_fqM.basicM <- function(t, y, pars, s) {
  with(pars$MYZpar[[s]], f*q)*y[pars$ix$MYZ[[s]]$M_ix]
}

#' @title Number of eggs laid by adult mosquitoes
#' @description Implements [F_eggs] for the basicM model.
#' @inheritParams F_eggs
#' @return a [numeric] vector of length `nPatches`
#' @export
F_eggs.basicM <- function(t, y, pars, s) {
  M <- y[pars$ix$MYZ[[s]]$M_ix]
  with(pars$MYZpar[[s]],{
    return(M*nu*eggsPerBatch)
  })
}

#' @title Derivatives for adult mosquitoes
#' @description Implements [dMYZdt] for the basicM model.
#' @inheritParams dMYZdt
#' @return a [numeric] vector
#' @export
dMYZdt.basicM <- function(t, y, pars, s) {
  Lambda = pars$Lambda[[s]]

  with(list_MYZvars(y, pars, s),{
    with(pars$MYZpar[[s]],{

      Mt <- Lambda + Omega %*% M
      Pt <- f*(M-P) + Omega %*% P

      return(c(Mt, Pt))
    })
  })
}

#' @title Setup MYZpar for the basicM model
#' @description Implements [setup_MYZpar] for the basicM model
#' @inheritParams setup_MYZpar
#' @return a [list] vector
#' @export
setup_MYZpar.basicM = function(MYZname, pars, s, MYZopts=list(), EIPname, calK){
  pars$MYZpar[[s]] = make_MYZpar_basicM(pars$nPatches, MYZopts, EIPname, calK)
  return(pars)
}


#' @title Make parameters for basicM adult mosquito model
#' @param nPatches is the number of patches, an integer
#' @param MYZopts a [list] of values that overwrites the defaults
#' @param EIPname a string: the class name for the EIP model
#' @param calK a mosquito dispersal matrix of dimensions `nPatches` by `nPatches`
#' @param p daily mosquito survival
#' @param sigma emigration rate
#' @param f feeding rate
#' @param q human blood fraction
#' @param eip the maximum number of cohorts in the EIP
#' @param nu oviposition rate, per mosquito
#' @param eggsPerBatch eggs laid per oviposition
#' @param p_mod a name to dispatch F_p
#' @param sigma_mod a name to dispatch F_sigma
#' @param f_mod a name to dispatch F_f
#' @param q_mod a name to dispatch F_q
#' @param nu_mod a name to dispatch F_nu
#' @return a [list]
#' @export
make_MYZpar_basicM = function(nPatches, MYZopts=list(), EIPname, calK,
                          p=11/12, sigma=1/8, f=0.3, q=0.95, eip=12,
                          nu=1, eggsPerBatch=60,
                          p_mod = "static",
                          sigma_mod = "static",
                          f_mod = "static",
                          q_mod = "static",
                          nu_mod = "static"){

  stopifnot(is.matrix(calK))
  stopifnot(dim(calK) == c(nPatches, nPatches))

  with(MYZopts,{
    MYZpar <- list()
    class(MYZpar) <- "basicM"

    MYZpar$nPatches <- nPatches
    if(nPatches == 1){
      sigma = 0
      calK = 1
    }

    MYZpar$p       <- checkIt(p, nPatches)
    MYZpar$sigma   <- checkIt(sigma, nPatches)
    MYZpar$f       <- checkIt(f, nPatches)
    MYZpar$q       <- checkIt(q, nPatches)
    MYZpar$nu      <- checkIt(nu, nPatches)
    MYZpar$eggsPerBatch <- eggsPerBatch

    # Store as baseline values
    MYZpar$p0      <- MYZpar$p
    MYZpar$sigma0  <- MYZpar$sigma
    MYZpar$f0      <- MYZpar$f
    MYZpar$q0      <- MYZpar$q
    MYZpar$nu0     <- MYZpar$nu

    MYZpar$p_par   <- list()
    class(MYZpar$p_par) <- p_mod
    MYZpar$f_par   <- list()
    class(MYZpar$f_par) <- f_mod
    MYZpar$q_par   <- list()
    class(MYZpar$q_par) <- q_mod
    MYZpar$sigma_par   <- list()
    class(MYZpar$sigma_par) <- sigma_mod
    MYZpar$nu_par   <- list()
    class(MYZpar$nu_par) <- nu_mod

    MYZpar$calK <- calK
    MYZpar$Omega <- make_Omega(0, MYZpar)

    return(MYZpar)
})}

#' @title Setup initial values for the basicM model
#' @description Implements [setup_MYZinits] for the basicM model
#' @inheritParams setup_MYZinits
#' @return a [list]
#' @export
setup_MYZinits.basicM = function(pars, s, MYZopts=list()){
  pars$MYZinits[[s]] = with(pars$MYZpar[[s]], make_MYZinits_basicM(nPatches, MYZopts))
  return(pars)
}

#' @title Make inits for basicM adult mosquito model
#' @param nPatches the number of patches in the model
#' @param MYZopts a [list] of values that overwrites the defaults
#' @param M0 total mosquito density at each patch
#' @param P0 total parous mosquito density at each patch
#' @return a [list]
#' @export
make_MYZinits_basicM = function(nPatches, MYZopts = list(),
                            M0=5, P0=1){
  with(MYZopts,{
    M = checkIt(M0, nPatches)
    P = checkIt(P0, nPatches)
    return(list(M=M, P=P))
  })
}


#' @title Add indices for adult mosquitoes to parameter list
#' @description Implements [make_indices_MYZ] for the basicM model.
#' @inheritParams make_indices_MYZ
#' @return a [list]
#' @importFrom utils tail
#' @export
make_indices_MYZ.basicM <- function(pars, s) {with(pars,{

  M_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(M_ix, 1)

  P_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(P_ix, 1)

  pars$max_ix = max_ix
  pars$ix$MYZ[[s]] = list(M_ix=M_ix, P_ix=P_ix)
  return(pars)
})}


#' @title Return the variables as a list
#' @description This method dispatches on the type of `pars$MYZpar[[s]]`
#' @inheritParams list_MYZvars
#' @return a [list]
#' @export
list_MYZvars.basicM <- function(y, pars, s){
  with(pars$ix$MYZ[[s]],
       return(list(
         M = y[M_ix],
         P = y[P_ix],
       )))
}

#' @title Make parameters for basicM adult mosquito model
#' @param pars a [list]
#' @param p daily mosquito survival
#' @param sigma emigration rate
#' @param f feeding rate
#' @param q human blood fraction
#' @param nu oviposition rate, per mosquito
#' @param eggsPerBatch eggs laid per oviposition
#' @param calK mosquito dispersal matrix of dimensions `nPatches` by `nPatches`
#' @return a [list]
#' @export
make_parameters_MYZ_basicM <- function(pars,  p, sigma, f, q, nu, eggsPerBatch,  calK) {
  stopifnot(is.numeric(p), is.numeric(sigma), is.numeric(f),
            is.numeric(q), is.numeric(nu), is.numeric(eggsPerBatch))

  MYZpar <- list()
  class(MYZpar) <- "basicM"

  nPatches <- pars$nPatches
  MYZpar$nPatches <- nPatches
  if(nPatches == 1){
    sigma = 0
    calK = 1
  }

  MYZpar$p      <- checkIt(p, pars$nPatches)
  MYZpar$sigma  <- checkIt(sigma, pars$nPatches)
  MYZpar$f      <- checkIt(f, pars$nPatches)
  MYZpar$q      <- checkIt(q, pars$nPatches)
  MYZpar$nu     <- checkIt(nu, pars$nPatches)
  MYZpar$eggsPerBatch <- eggsPerBatch

  # Store as baseline values
  MYZpar$p0      <- MYZpar$p
  MYZpar$sigma0  <- MYZpar$sigma
  MYZpar$f0      <- MYZpar$f
  MYZpar$q0      <- MYZpar$q
  MYZpar$nu0     <- MYZpar$nu


  MYZpar$nPatches <- pars$nPatches

  MYZpar$calK <- calK
  Omega   <- make_Omega(0, MYZpar)

  pars$MYZpar = list()
  pars$MYZpar[[1]] = MYZpar

  return(pars)
}

#' @title Make inits for basicM adult mosquito model
#' @param pars a [list]
#' @param M0 total mosquito density at each patch
#' @param P0 total parous mosquito density at each patch
#' @return a [list]
#' @export
make_inits_MYZ_basicM <- function(pars, M0, P0) {
  pars$MYZinits = list()
  pars$MYZinits[[1]] = list(M=M0, P=P0)
  return(pars)
}

#' @title Parse the output of deSolve and return variables for the basicM model
#' @description Implements [parse_dts_out_MYZ] for the basicM model
#' @inheritParams parse_dts_out_MYZ
#' @return a [list]
#' @export
parse_dts_out_MYZ.basicM <- function(dts_out, pars, s) {with(pars$ix$MYZ[[s]],{
  time = dts_out[,1]
  M = dts_out[,M_ix+1]
  P = dts_out[,P_ix+1]
  parous = P/M
  return(list(time=time, M=M, P=P, parous=parous))
})}

#' @title Return initial values as a vector
#' @description Implements [get_inits_MYZ] for the basicM model.
#' @inheritParams get_inits_MYZ
#' @return [numeric]
#' @export
get_inits_MYZ.basicM <- function(pars, s) {with(pars$MYZinits[[s]],{
  c(M, P)
})}

#' @title Make inits for basicM adult mosquito model
#' @inheritParams update_inits_MYZ
#' @return a [list]
#' @export
update_inits_MYZ.basicM <- function(pars, y0, s) {with(pars$ix$MYZ[[s]],{
  M = y0[M_ix]
  P = y0[P_ix]
  pars$MYZinits[[s]] = make_MYZinits_basicM(pars$nPatches,  list(), M0=M, P=P0)
  return(pars)
})}


