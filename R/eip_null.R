
#' @title Modify parameters due to exogenous forcing by all kinds of control
#' @description Implements [EIP] for the null model
#' @inheritParams EIP
#' @return [numeric]
#' @export
EIP.null <- function(t, MYZpar){numeric(0)}

#' @title Set up a null model for the EIP
#' @inheritParams setup_EIP
#' @return [list]
#' @export
setup_EIP.null<- function(EIPname, MYZpar, MYZopts=list()){
  setup_eip_null(MYZopts, MYZpar)
}

#' @title Set up a null model for the EIP
#' @param MYZopts a [list]
#' @param MYZpar the MYZ parameters
#' @return [list]
#' @export
setup_eip_null = function(MYZopts=list(), MYZpar){
  EIPmod <- list()
  class(EIPmod) <- 'null'
  MYZpar$EIPmod <- EIPmod
  return(MYZpar)
}
