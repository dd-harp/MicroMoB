
#' @title Modify parameters due to exogenous forcing by all kinds of control
#' @description Implements [EIP] for the fixed model (the EIP is constant)
#' @inheritParams EIP
#' @return the EIP maturation vector G [numeric]
#' @export
EIP.fixed <- function(t, MYZpar){
  G <- rep(0, MYZpar$max_eip)
  G[MYZpar$eip] <- 1
  return(G)
}

#' @title Set up a fixed model for the EIP
#' @inheritParams setup_EIP
#' @return [list]
#' @export
setup_EIP.fixed <- function(EIPname, MYZpar, MYZopts=list()){
  setup_eip_fixed(MYZopts, MYZpar)
}

#' @title Set up a fixed model for the EIP
#' @param MYZopts a [list]
#' @param MYZpar the MYZ parameters
#' @return [list]
#' @export
setup_eip_fixed = function(MYZopts=list(), MYZpar){
  with(MYZpar,
    with(MYZopts,{
      EIPmod <- list()
      class(EIPmod) <- 'fixed'
      MYZpar$EIPmod <- EIPmod
      MYZpar$eip <- eip
      MYZpar$max_eip <- eip
      return(MYZpar)
}))}
