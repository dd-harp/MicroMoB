# specialized methods for the adult mosquito RM model

#' @title Reset bloodfeeding and mortality rates to baseline
#' @description Implements [MBionomics] for the RM model
#' @inheritParams MBionomics
#' @return the model as a [list]
#' @export
MBionomics.RM <- function(t, y, pars, s) {
  with(pars$MYZpar[[s]],{
    pars$MYZpar[[s]]$f <- f0
    pars$MYZpar[[s]]$q <- q0
    pars$MYZpar[[s]]$p <- p0
    pars$MYZpar[[s]]$sigma <- sigma0
    pars$MYZpar[[s]]$nu <- nu0
    #pars$MYZpar[[s]]$eip <- EIP(t, EIPmod)

    return(pars)
})}

#' @title The net blood feeding rate of the infective mosquito population in a patch
#' @description Implements [F_fqZ] for the RM model.
#' @inheritParams F_fqZ
#' @return a [numeric] vector of length `nPatches`
#' @export
F_fqZ.RM <- function(t, y, pars, s) {
  with(pars$MYZpar[[s]], f*q)*y[pars$ix$MYZ[[s]]$Z_ix]
}

#' @title The net blood feeding rate of the infective mosquito population in a patch
#' @description Implements [F_fqM] for the RM model.
#' @inheritParams F_fqM
#' @return a [numeric] vector of length `nPatches`
#' @export
F_fqM.RM <- function(t, y, pars, s) {
  with(pars$MYZpar[[s]], f*q)*y[pars$ix$MYZ[[s]]$Z_ix]
}

#' @title Number of eggs laid by adult mosquitoes
#' @description Implements [F_eggs] for the RM model.
#' @inheritParams F_eggs
#' @return a [numeric] vector of length `nPatches`
#' @export
F_eggs.RM <- function(t, y, pars, s) {
  M <- y[pars$ix$MYZ[[s]]$M_ix]
  with(pars$MYZpar[[s]],{
    return(M*nu*eggsPerBatch)
  })
}

#' @title Derivatives for adult mosquitoes
#' @description Implements [dMYZdt] for the RM ODE model.
#' @inheritParams dMYZdt
#' @return a [numeric] vector
#' @export
dMYZdt.RM <- function(t, y, pars, s) {
  Lambda = pars$Lambda[[s]]
  kappa = pars$kappa[[s]]

  with(list_MYZvars(y, pars, s),{
    with(pars$MYZpar[[s]],{
        Omega <- make_Omega(p, sigma, calK, nPatches)
        #G <- EIP(t, EIPmod)

        Mt <- Lambda + Omega %*% M
        Pt <- f*(M-P) + Omega %*% P
        Ut <- Lambda + Omega %*% (exp(-f*q*kappa)*U)
        Yt <- Omega %*% (Y %*% diag(1-G))
        Zt <- Omega %*% (Y%*%G)  + (Omega %*% Z)

        eip_day_ix = (t %% max_eip) + 1
        eip_yday_ix = ((t-1) %% max_eip) + 1

        Yt[,eip_yday_ix]  <- Yt[,eip_yday_ix] + Yt[,eip_day_ix]
        Yt[,eip_day_ix] <- Omega %*% ((1-exp(-f*q*kappa))*U)

        return(c(Mt, Pt, Ut, Yt, Zt))
      })
  })
}

#' @title Setup MYZpar for the RM model
#' @description Implements [setup_MYZpar] for the RM model
#' @inheritParams setup_MYZpar
#' @return a [list] vector
#' @export
setup_MYZpar.RM = function(MYZname, pars, s, MYZopts=list(), EIPmod, calK){
  pars$MYZpar[[s]] = make_MYZpar_RM(pars$nPatches, MYZopts, EIPmod, calK)
  return(pars)
}

#' @title Make parameters for RM ODE adult mosquito model
#' @param nPatches is the number of patches, an integer
#' @param MYZopts a [list] of values that overwrites the defaults
#' @param EIPmod a [list] that defines the EIP model
#' @param calK a mosquito dispersal matrix of dimensions `nPatches` by `nPatches`
#' @param p daily mosquito survival
#' @param sigma emigration rate
#' @param f feeding rate
#' @param q human blood fraction
#' @param eip the maximum number of cohorts in the EIP
#' @param nu oviposition rate, per mosquito
#' @param eggsPerBatch eggs laid per oviposition
#' @return a [list]
#' @export
make_MYZpar_RM = function(nPatches, MYZopts=list(), EIPmod, calK,
                          p=11/12, sigma=1/8, f=0.3, q=0.95, eip=12,
                          nu=1, eggsPerBatch=60){

  stopifnot(is.matrix(calK))
  stopifnot(dim(calK) == c(nPatches, nPatches))

  with(MYZopts,{
    MYZpar <- list()
    class(MYZpar) <- "RM"

    MYZpar$nPatches <- nPatches

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

    # The EIP model and the eip
    MYZpar$eip <- eip
    MYZpar$max_eip <- eip
    MYZpar$G = rep(0, eip)
    MYZpar$G[eip] = 1

    #MYZpar$EIPmod <- EIPmod
    #MYZpar$eip <- EIP(0, EIPmod)

    MYZpar$calK <- calK

    MYZpar$Omega <- make_Omega(p, sigma, calK, nPatches)

    return(MYZpar)
})}

#' @title Setup initial values for the RM model
#' @description Implements [setup_MYZinits] for the RM model
#' @inheritParams setup_MYZinits
#' @return a [list]
#' @export
setup_MYZinits.RM = function(pars, s, MYZopts=list()){
  pars$MYZinits[[s]] = with(pars$MYZpar[[s]], make_MYZinits_RM(nPatches, max_eip, MYZopts))
  return(pars)
}

#' @title Make inits for RM adult mosquito model
#' @param nPatches the number of patches in the model
#' @param max_eip the maximum number of EIP cohorts, an [integer]
#' @param MYZopts a [list] of values that overwrites the defaults
#' @param M0 total mosquito density at each patch
#' @param P0 total parous mosquito density at each patch
#' @param U0 total uninfected mosquito density at each patch
#' @param Y0 infected mosquito density at each patch
#' @param Z0 infectious mosquito density at each patch
#' @return a [list]
#' @export
make_MYZinits_RM = function(nPatches, max_eip, MYZopts = list(),
                                M0=5, P0=1, U0=0, Y0=1, Z0=1){
  with(MYZopts,{
    M = checkIt(M0, nPatches)
    P = checkIt(P0, nPatches)
    U = checkIt(U0, nPatches)
    Y = checkIt(Y0, nPatches*max_eip)
    Z = checkIt(Z0, nPatches)
    return(list(M=M, P=P, U=U, Y=Y, Z=Z))
  })
}


#' @title Add indices for adult mosquitoes to parameter list
#' @description Implements [make_indices_MYZ] for the RM model.
#' @inheritParams make_indices_MYZ
#' @return a [list]
#' @importFrom utils tail
#' @export
make_indices_MYZ.RM <- function(pars, s) {with(pars,{

  M_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(M_ix, 1)

  P_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(P_ix, 1)

  U_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(U_ix, 1)

  Y_ix <- seq(from = max_ix+1, length.out=nPatches*MYZpar[[s]]$max_eip)
  max_ix <- tail(Y_ix, 1)

  Z_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(Z_ix, 1)

  pars$max_ix = max_ix
  pars$ix$MYZ[[s]] = list(M_ix=M_ix, P_ix=P_ix, U_ix=U_ix, Y_ix=Y_ix, Z_ix=Z_ix)
  return(pars)
})}


#' @title Return the variables as a list
#' @description This method dispatches on the type of `pars$MYZpar[[s]]`
#' @inheritParams list_MYZvars
#' @return a [list]
#' @export
list_MYZvars.RM <- function(y, pars, s){
  with(pars$ix$MYZ[[s]],
       return(list(
         M = y[M_ix],
         P = y[P_ix],
         U = y[U_ix],
         Y = y[Y_ix],
         Z = y[Z_ix]
  )))
}

#' @title Make parameters for RM ODE adult mosquito model
#' @param pars a [list]
#' @param p daily mosquito survival
#' @param sigma emigration rate
#' @param f feeding rate
#' @param q human blood fraction
#' @param nu oviposition rate, per mosquito
#' @param eggsPerBatch eggs laid per oviposition
#' @param eip maximum length for extrinsic incubation period
#' @param calK mosquito dispersal matrix of dimensions `nPatches` by `nPatches`
#' @return a [list]
#' @export
make_parameters_MYZ_RM <- function(pars, p, sigma, f, q, nu, eggsPerBatch, eip, calK) {
  stopifnot(is.numeric(p), is.numeric(sigma), is.numeric(f),
            is.numeric(q), is.numeric(nu), is.numeric(eggsPerBatch))

  MYZpar <- list()
  class(MYZpar) <- "RM"

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

  Omega   <- expm::expm(-make_Omega(p, sigma, calK, pars$nPatches))
  MYZpar$eip <- eip
  MYZpar$max_eip <- eip
  MYZpar$G <- rep(0, eip)
  MYZpar$G[eip] <- 1
  #MYZpar$EIPmod <- setup_eip_static(eip=eip)
  #MYZpar$eip <- eip
  MYZpar$calK <- calK
  MYZpar$nPatches <- pars$nPatches

  pars$MYZpar = list()
  pars$MYZpar[[1]] = MYZpar

  return(pars)
}

#' @title Make inits for RM adult mosquito model
#' @param pars a [list]
#' @param M0 total mosquito density at each patch
#' @param P0 total parous mosquito density at each patch
#' @param U0 total uninfected mosquito density at each patch
#' @param Y0 infected mosquito density at each patch
#' @param Z0 infectious mosquito density at each patch
#' @return a [list]
#' @export
make_inits_MYZ_RM <- function(pars, M0, P0, U0, Y0, Z0) {
  pars$MYZinits = list()
  pars$MYZinits[[1]] = list(M=M0, P=P0, U=U0, Y=Y0, Z=Z0)
  return(pars)
}

#' @title Parse the output of deSolve and return variables for the RM model
#' @description Implements [parse_dts_out_MYZ] for the RM model
#' @inheritParams parse_dts_out_MYZ
#' @return a [list]
#' @export
parse_dts_out_MYZ.RM <- function(dts_out, pars, s) {with(pars$ix$MYZ[[s]],{
  time = dts_out[,1]
  M = dts_out[,M_ix+1]
  P = dts_out[,P_ix+1]
  U = dts_out[,U_ix+1]
  Y = rowSums(dts_out[,Y_ix+1])
  Z = dts_out[,Z_ix+1]
  y = Y/M
  z = Z/M
  parous = P/M
  return(list(time=time, M=M, P=P, U=U, Y=Y, Z=Z, y=y, z=z, parous=parous))
})}

#' @title Return initial values as a vector
#' @description Implements [get_inits_MYZ] for the RM model.
#' @inheritParams get_inits_MYZ
#' @return [numeric]
#' @export
get_inits_MYZ.RM <- function(pars, s) {with(pars$MYZinits[[s]],{
  c(M, P, U, Y, Z)
})}

#' @title Make inits for RM adult mosquito model
#' @inheritParams update_inits_MYZ
#' @return a [list]
#' @export
update_inits_MYZ.RM <- function(pars, y0, s) {with(pars$ix$MYZ[[s]],{
  M = y0[M_ix]
  P = y0[P_ix]
  U = y0[U_ix]
  Y = y0[Y_ix]
  Z = y0[Z_ix]
  pars$MYZinits[[s]] = make_MYZinits_RM(pars$nPatches, max_eip, list(), M0=M, P=P0, U=U0, Y0=Y, Z0=Z)
  return(pars)
})}



