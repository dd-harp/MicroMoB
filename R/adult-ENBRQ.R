# specialized methods for the adult mosquito ENBRQ model

#' @title Reset bloodfeeding and mortality rates to baseline
#' @description Implements [MBionomics] for the ENBRQ model
#' @inheritParams MBionomics
#' @return the model as a [list]
#' @export
MBionomics.ENBRQ <- function(t, y, pars, s) {
  pars$MYZpar[[s]]$f     <- F_f(t, pars$MYZpar[[s]])
  pars$MYZpar[[s]]$q     <- F_q(t, pars$MYZpar[[s]])
  pars$MYZpar[[s]]$p     <- F_p(t, pars$MYZpar[[s]])
  pars$MYZpar[[s]]$sigma <- F_sigma(t, pars$MYZpar[[s]])
  pars$MYZpar[[s]]$nu    <- F_nu(t, pars$MYZpar[[s]])
  return(pars)
}

#' @title The net blood feeding rate of the infective mosquito population in a patch
#' @description Implements [F_fqZ] for the ENBRQ model.
#' @inheritParams F_fqZ
#' @return a [numeric] vector of length `nPatches`
#' @export
F_fqZ.ENBRQ <- function(t, y, pars, s) {
  numeric(0)
}

#' @title The net blood feeding rate of the infective mosquito population in a patch
#' @description Implements [F_fqM] for the ENBRQ model.
#' @inheritParams F_fqM
#' @return a [numeric] vector of length `nPatches`
#' @export
F_fqM.ENBRQ <- function(t, y, pars, s) {
  y[pars$ix$MYZ[[s]]$fqM_ix]
}

#' @title Number of eggs laid by adult mosquitoes
#' @description Implements [F_eggs] for the ENBRQ model.
#' @inheritParams F_eggs
#' @return a [numeric] vector of length `nPatches`
#' @export
F_eggs.ENBRQ <- function(t, y, pars, s) {
  y[pars$ix$MYZ[[s]]$eggs_ix]
}

#' @title Derivatives for adult mosquitoes
#' @description Implements [dMYZdt] for the ENBRQ model.
#' @inheritParams dMYZdt
#' @return a [numeric] vector
#' @export
dMYZdt.ENBRQ <- function(t, y, pars, s) {
  Lambda = pars$Lambda[[s]]
  with(list_MYZvars(y, pars, s),{
    with(pars$MYZpar[[s]],{

      fqM = 0*B
      eggs = 0*Q
      for(d in 1:D){
        dix = d%%nR1 + 1

        Et <- Lambda/D
        Nt <- K_bq %*% (p[,d]*E)        + K_bb %*% ((1-f[,d])*p*N)
        Bt <- K_bb %*% ((1-f[,d])*p[,d]*B)  + K_bq %*% (nu[,d]*p[,d]*Q)
        Qt <- K_qr %*% (r2q*p[,d]*R2t) + K_bq %*% (p[,d]*(1-nu[,d])*Q)

        R1t <- p[,d]*R1t
        R2t <- p[,d]*(1-r2q)*R1t[,dix]
        R1t[,dix] <- K_rb %*% (f[,d]*p[,d]*(N+B))

        fqM = fqM + f[,d]*q[,d]*(N+B)
        eggs = eggs + nu[,d]*eggsPerBatch*Q

        E=Et; N=Nt; B=Bt; Q=Qt; R1=R1t; R2=R2t
      }
      return(c(Et, Nt, Bt, Qt, R1t, R2t, fqM, eggs))
    })
  })
}

#' @title Setup MYZpar for the ENBRQ model
#' @description Implements [setup_MYZpar] for the ENBRQ model
#' @inheritParams setup_MYZpar
#' @return a [list] vector
#' @export
setup_MYZpar.ENBRQ = function(MYZname, pars, s, MYZopts=list(), EIPname, calK){
  pars$MYZpar[[s]] = make_MYZpar_ENBRQ(pars$nPatches, MYZopts, calK)
  return(pars)
}


#' @title Make parameters for ENBRQ adult mosquito model
#' @param nPatches is the number of patches, an integer
#' @param MYZopts a [list] of values that overwrites the defaults
#' @param calK a mosquito dispersal matrix of dimensions `nPatches` by `nPatches`
#' @param D number of time steps per day
#' @param nR1 number of time steps in R1
#' @param p daily mosquito survival
#' @param sigma emigration rate
#' @param f feeding rate
#' @param q human blood fraction
#' @param r2q proportion of R2 transitioning to Q
#' @param nu oviposition rate, per mosquito
#' @param eggsPerBatch eggs laid per oviposition
#' @param p_mod a name to dispatch F_p
#' @param sigma_mod a name to dispatch F_sigma
#' @param f_mod a name to dispatch F_f
#' @param q_mod a name to dispatch F_q
#' @param nu_mod a name to dispatch F_nu
#' @return a [list]
#' @export
make_MYZpar_ENBRQ = function(nPatches, MYZopts=list(),  calK, D=4, nR1=3,
                          p=11/12,
                          sigma=1/8,
                          f=0.3,
                          q=0.95,
                          r2q=0.5,
                          nu=1,
                          eggsPerBatch=60,
                          p_mod = "dddn",
                          sigma_mod = "dddn",
                          f_mod = "dddn",
                          q_mod = "dddn",
                          nu_mod = "dddn"
){

  stopifnot(is.matrix(calK))
  stopifnot(dim(calK) == c(nPatches, nPatches))

  with(MYZopts,{
    MYZpar <- list()
    class(MYZpar) <- "ENBRQ"

    MYZpar$nPatches <- nPatches
    if(nPatches == 1){
      sigma = 0
      calK = 1
    }

    MYZpar$D       <- checkIt(D, 1)
    MYZpar$nR1    <- checkIt(nR1, 1)
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
    class(MYZpar$p_par) <- "dddn"
    MYZpar$f_par   <- list()
    class(MYZpar$f_par) <- "dddn"
    MYZpar$q_par   <- list()
    class(MYZpar$q_par) <-  "dddn"
    MYZpar$sigma_par   <- list()
    class(MYZpar$sigma_par) <- "dddn"
    MYZpar$nu_par   <- list()
    class(MYZpar$nu_par) <- "dddn"

    MYZpar$K_bb <- calK
    MYZpar$K_bq <- calK
    MYZpar$K_qb <- calK
    MYZpar$K_qq <- calK
    MYZpar$K_qr <- calK
    MYZpar$K_rb <- calK

    #MYZpar$K_bb <- with(calK, K_bb)
    #MYZpar$K_bq <- with(calK, K_bq)
    #MYZpar$K_qb <- with(calK, K_qb)
    #MYZpar$K_qq <- with(calK, K_qq)
    #MYZpar$K_qr <- with(calK, K_qr)
    #MYZpar$K_rb <- with(calK, K_rb)

    return(MYZpar)
  })}

#' @title Setup initial values for the ENBRQ model
#' @description Implements [setup_MYZinits] for the ENBRQ model
#' @inheritParams setup_MYZinits
#' @return a [list]
#' @export
setup_MYZinits.ENBRQ = function(pars, s, MYZopts=list()){
  pars$MYZinits[[s]] = with(pars$MYZpar[[s]], make_MYZinits_ENBRQ(nPatches, nR1, MYZopts))
  return(pars)
}

#' @title Make inits for ENBRQ adult mosquito model
#' @param nPatches the number of patches in the model
#' @param nR1 number of time steps in R1
#' @param MYZopts a [list] of values that overwrites the defaults
#' @param E0 total null-parous, blood feeding mosquito density
#' @param N0 total null-parous, blood feeding mosquito density
#' @param B0 total parous, blood feeding mosquito density at each patch
#' @param Q0 total egg laying mosquito density at each patch
#' @param R10 resting, stage 1, mosquito density at each patch
#' @param R20 resting, stage 2, mosquito density at each patch
#' @return a [list]
#' @export
make_MYZinits_ENBRQ = function(nPatches, nR1, MYZopts = list(),
                            E0=5, N0=5, B0=5, Q0=1, R10=0, R20=0){
  with(MYZopts,{
    E = checkIt(E0, nPatches)
    N = checkIt(N0, nPatches)
    B = checkIt(B0, nPatches)
    Q = checkIt(Q0, nPatches)
    R1 = checkIt(R10, nPatches*nR1)
    R2 = checkIt(R20, nPatches)
    return(list(E=E, N=N, B=B, Q=Q, R1=R1, R2=R2))
  })
}


#' @title Add indices for adult mosquitoes to parameter list
#' @description Implements [make_indices_MYZ] for the ENBRQ model.
#' @inheritParams make_indices_MYZ
#' @return a [list]
#' @importFrom utils tail
#' @export
make_indices_MYZ.ENBRQ <- function(pars, s) {with(pars,{

  E_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(E_ix, 1)

  N_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(N_ix, 1)

  B_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(B_ix, 1)

  Q_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(Q_ix, 1)

  R1_ix <- seq(from = max_ix+1, length.out=nPatches*pars$MYZpar[[s]]$nR1)
  max_ix <- tail(R1_ix, 1)

  R2_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(R2_ix, 1)

  fqM_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(fqM_ix, 1)

  eggs_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(eggs_ix, 1)

  pars$max_ix = max_ix
  pars$ix$MYZ[[s]] = list(E_ix=E_ix, N_ix=N_ix, B_ix=B_ix, Q_ix=Q_ix, R1_ix = R1_ix, R2_ix=R2_ix, fqM_ix=fqM_ix, eggs_ix=eggs_ix)
  return(pars)
})}


#' @title Return the variables as a list
#' @description This method dispatches on the type of `pars$MYZpar[[s]]`
#' @inheritParams list_MYZvars
#' @return a [list]
#' @export
list_MYZvars.ENBRQ <- function(y, pars, s){
  with(pars$ix$MYZ[[s]],
       return(list(
         E = y[E_ix],
         N = y[N_ix],
         B = y[B_ix],
         Q = y[Q_ix],
         R1 = with(pars$MYZpar[[s]], matrix(y[R1_ix], nR1, nPatches)),
         R2 = y[R2_ix]
       )))
}

#' @title Make parameters for ENBRQ adult mosquito model
#' @param pars a [list]
#' @param D time steps per day
#' @param nR1 time steps per day
#' @param p daily mosquito survival
#' @param sigma emigration rate
#' @param f feeding rate
#' @param q human blood fraction
#' @param r2q proportion of R2 transitioning to Q
#' @param nu oviposition rate, per mosquito
#' @param eggsPerBatch eggs laid per oviposition
#' @param calK mosquito dispersal matrix of dimensions `nPatches` by `nPatches`
#' @return a [list]
#' @export
make_parameters_MYZ_ENBRQ <- function(pars, D, nR1, p, sigma, f, q, r2q, nu, eggsPerBatch, calK) {
  stopifnot(is.numeric(p), is.numeric(sigma), is.numeric(f),
            is.numeric(q), is.numeric(nu), is.numeric(eggsPerBatch))

  MYZpar <- list()
  class(MYZpar) <- "ENBRQ"

  nPatches <- pars$nPatches
  MYZpar$nPatches <- nPatches
  if(nPatches == 1){
    sigma = 0
    calK = 1
  }

  MYZpar$D      <- checkIt(D, 1)
  MYZpar$nR1    <- checkIt(nR1, 1)
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

  MYZpar$p_par   <- list()
  class(MYZpar$p_par) <- "dddn"
  MYZpar$f_par   <- list()
  class(MYZpar$f_par) <- "dddn"
  MYZpar$q_par   <- list()
  class(MYZpar$q_par) <- "dddn"
  MYZpar$sigma_par   <- list()
  class(MYZpar$sigma_par) <- "dddn"
  MYZpar$nu_par   <- list()
  class(MYZpar$nu_par) <- "dddn"

  MYZpar$nPatches <- pars$nPatches

  MYZpar$K_bb <- calK
  MYZpar$K_bq <- calK
  MYZpar$K_qb <- calK
  MYZpar$K_qq <- calK
  MYZpar$K_qr <- calK
  MYZpar$K_rb <- calK

  #MYZpar$K_bb <- with(calK, K_bb)
  #MYZpar$K_bq <- with(calK, K_bq)
  #MYZpar$K_qb <- with(calK, K_qb)
  #MYZpar$K_qq <- with(calK, K_qq)
  #MYZpar$K_qr <- with(calK, K_qr)
  #MYZpar$K_rb <- with(calK, K_rb)

  pars$MYZpar = list()
  pars$MYZpar[[1]] = MYZpar

  return(pars)
}

#' @title Make inits for ENBRQ adult mosquito model
#' @param pars a [list]
#' @param E0 total mosquito density at each patch
#' @param N0 total mosquito density at each patch
#' @param B0 blood feeding mosquito density at each patch
#' @param Q0 egg laying mosquito density at each patch
#' @param R10 total stage 1 resting mosquitoes in each patch
#' @param R20 total stage 2 resting mosquitoes in each patch
#' @return a [list]
#' @export
make_inits_MYZ_ENBRQ <- function(pars, E0, N0, B0, Q0, R10, R20) {
  pars$MYZinits = list()
  pars$MYZinits[[1]] = list(E=E0, N=N0, B=B0, Q=Q0, R1=R10, R2=R20)
  return(pars)
}

#' @title Parse the output of deSolve and return variables for the ENBRQ model
#' @description Implements [parse_dts_out_MYZ] for the ENBRQ model
#' @inheritParams parse_dts_out_MYZ
#' @return a [list]
#' @export
parse_dts_out_MYZ.ENBRQ <- function(dts_out, pars, s) {with(pars$ix$MYZ[[s]],{
  time = dts_out[,1]
  E = dts_out[,E_ix+1]
  N = dts_out[,N_ix+1]
  B = dts_out[,B_ix+1]
  Q = dts_out[,Q_ix+1]
  R1 = colSums(with(pars$MYZpar[[s]], matrix(dts_out[,R1_ix+1], nR1, nPatches)))
  R2 = dts_out[,R2_ix+1]
  fqM = dts_out[,fqM_ix+1]
  eggs = dts_out[,eggs_ix+1]
  return(list(time=time, N=N, B=B, Q=Q, R1=R1, R2=R2, fqM=fqM, eggs=eggs))
})}

#' @title Return initial values as a vector
#' @description Implements [get_inits_MYZ] for the ENBRQ model.
#' @inheritParams get_inits_MYZ
#' @return [numeric]
#' @export
get_inits_MYZ.ENBRQ <- function(pars, s) {with(pars$MYZinits[[s]],{
  c(E, N, B, Q, R1, R2)
})}

#' @title Make inits for ENBRQ adult mosquito model
#' @inheritParams update_inits_MYZ
#' @return a [list]
#' @export
update_inits_MYZ.ENBRQ <- function(pars, y0, s) {with(pars$ix$MYZ[[s]],{
  E = y0[E_ix]
  N = y0[N_ix]
  B = y0[B_ix]
  Q = y0[Q_ix]
  R1 = y0[R1_ix]
  R2 = y0[R2_ix]
  pars$MYZinits[[s]] = make_MYZinits_ENBRQ(pars$nPatches, pars$MYZpars[[s]]$nR1, list(), E0=E, N0=N, B0=B, Q0=Q, R10=R1, R20=R2)
  return(pars)
})}



