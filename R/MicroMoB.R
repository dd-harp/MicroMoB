#' MicroMoB: Microsimulation for mosquito-borne pathogens
#'
#' @description
#' Discrete time simulation of mosquito-borne pathogen transmission
#'
#' @docType package
#' @name MicroMoB
#'
#' @importFrom utils getFromNamespace
#' @importFrom httpuv startServer stopServer
#'
#' @useDynLib MicroMoB, .registration = TRUE
"_PACKAGE"

findPort <- getFromNamespace("findPort", "plumber")
