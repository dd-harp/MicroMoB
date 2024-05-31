
#' @title Compute the blood feeding rate, f
#' @description This method dispatches on the type of `MYZpar$f_par`. It should
#' set the values of the bionomic parameters to baseline values.
#' @param t current simulation time
#' @param MYZpar a [list]
#' @return a [numeric] vector of length `nPatches`
#' @export
F_f <- function(t, MYZpar) {
  UseMethod("F_f", MYZpar$f_par)
}

#' @title Static model for the blood feeding rate
#' @description Implements [F_f] for a static model
#' @inheritParams F_f
#' @return a [numeric] vector of length `nPatches`
#' @export
F_f.static <- function(t, MYZpar){
  MYZpar$f0
}

#' @title Dawn, day, dusk, night model for the blood feeding rate
#' @description Implements [F_f] for a dddn model
#' @inheritParams F_f
#' @return a [numeric] vector of length `nPatches`
#' @export
F_f.dddn <- function(t, MYZpar){
  with(MYZpar, t(matrix(f0, 4, nPatches)))
}

#' @title Type 2 functional response for the blood feeding rate
#' @description Implements [F_f] for a static model
#' @inheritParams F_f
#' @return a [numeric] vector of length `nPatches`
#' @export
F_f.type2 <- function(t, MYZpar){
  B = with(MYZpar, W + Visitors + O)
  with(MYZpar$f_par,{
    return(fx*sf*B/(1+sf*B))
  })
}

#' @title Compute the human blood fraction
#' @description This method dispatches on the type of `MYZpar$q_par`. It should
#' set the values of the bionomic parameters to baseline values.
#' @param t current simulation time
#' @param MYZpar a [list]
#' @return a [numeric] vector of length `nPatches`
#' @export
F_q <- function(t, MYZpar) {
  UseMethod("F_q", MYZpar$q_par)
}

#' @title Static model for human blood fraction
#' @description Implements [F_q] for a static model
#' @inheritParams F_q
#' @return a [numeric] vector of length `nPatches`
#' @export
F_q.static <- function(t, MYZpar){
  MYZpar$q0
}

#' @title Dawn, day, dusk, night model for the human fraction
#' @description Implements [F_q] for a dddn model
#' @inheritParams F_q
#' @return a [numeric] vector of length `nPatches`
#' @export
F_q.dddn <- function(t, MYZpar){
  with(MYZpar, t(matrix(q0, 4, nPatches)))
}

#' @title Static model for human blood fraction
#' @description Implements [F_q] for a static model
#' @inheritParams F_q
#' @return a [numeric] vector of length `nPatches`
#' @export
F_q.dynamic <- function(t, MYZpar){
  with(MYZpar,{
    return((W+Visitors)/(W + Visitors + O))
  })
}

#' @title Compute mosguito survival
#' @description This method dispatches on the type of `MYZpar$p_par`. It should
#' set the values of g to (possibly changing) baseline values.
#' @param t current simulation time
#' @param MYZpar a [list]
#' @return a [numeric] vector of length `nPatches`
#' @export
F_p <- function(t, MYZpar) {
  UseMethod("F_p", MYZpar$p_par)
}

#' @title Static model for mosquito survival
#' @description Implements [F_p] for a static model
#' @inheritParams F_p
#' @return a [numeric] vector of length `nPatches`
#' @export
F_p.static <- function(t, MYZpar){
  MYZpar$p0
}

#' @title Dawn, day, dusk, night model for the human fraction
#' @description Implements [F_p] for a dddn model
#' @inheritParams F_p
#' @return a [numeric] vector of length `nPatches`
#' @export
F_p.dddn <- function(t, MYZpar){
  with(MYZpar, t(matrix(p0, 4, nPatches)))
}

#' @title Compute mosquito emigration rates
#' @description This method dispatches on the type of `MYZpar$sigma_par`. It should
#' set the values of sigma to (possibly changing) baseline value(s).
#' @param t current simulation time
#' @param MYZpar a [list]
#' @return a [numeric] vector of length `nPatches`
#' @export
F_sigma <- function(t, MYZpar) {
  UseMethod("F_sigma", MYZpar$sigma_par)
}

#' @title Static model for mosquito emigration
#' @description Implements [F_sigma] for a static model
#' @inheritParams F_sigma
#' @return a [numeric] vector of length `nPatches`
#' @export
F_sigma.static <- function(t, MYZpar){
  MYZpar$sigma0
}

#' @title Dawn, day, dusk, night model for the human fraction
#' @description Implements [F_sigma] for a dddn model
#' @inheritParams F_sigma
#' @return a [numeric] vector of length `nPatches`
#' @export
F_sigma.dddn <- function(t, MYZpar){
  with(MYZpar, t(matrix(sigma0, 4, nPatches)))
}

#' @title Compute the egg laying rate
#' @description This method dispatches on the type of `MYZpar$nu_par`. It should
#' set the values of nu to (possibly changing) baseline value(s).
#' @param t current simulation time
#' @param MYZpar a [list]
#' @return a [numeric] vector of length `nPatches`
#' @export
F_nu <- function(t, MYZpar) {
  UseMethod("F_nu", MYZpar$nu_par)
}

#' @title Static model for the egg laying rate
#' @description Implements [F_nu] for a static model
#' @inheritParams F_nu
#' @return a [numeric] vector of length `nPatches`
#' @export
F_nu.static <- function(t, MYZpar){
  MYZpar$nu0
}

#' @title Dawn, day, dusk, night model for the human fraction
#' @description Implements [F_nu] for a dddn model
#' @inheritParams F_nu
#' @return a [numeric] vector of length `nPatches`
#' @export
F_nu.dddn <- function(t, MYZpar){
  with(MYZpar, t(matrix(nu0, 4, nPatches)))
}

#' @title Type 2 functional response for the blood feeding rate
#' @description Implements [F_nu] for a static model
#' @inheritParams F_nu
#' @return a [numeric] vector of length `nPatches`
#' @export
F_nu.type2 <- function(t, MYZpar){
  with(MYZpar$nu_par,{
    return(nux*snu*Q/(1+snu*Q))
  })
}
