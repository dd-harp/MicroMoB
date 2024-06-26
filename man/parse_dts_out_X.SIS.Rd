% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/human-SIS.R
\name{parse_dts_out_X.SIS}
\alias{parse_dts_out_X.SIS}
\title{Parse the output of deSolve and return variables for the SIS model}
\usage{
\method{parse_dts_out_X}{SIS}(dts_out, pars, i)
}
\arguments{
\item{dts_out}{a \link{matrix} of outputs from deSolve}

\item{pars}{a \link{list} that defines a model}

\item{i}{the host species index}
}
\value{
none
}
\description{
Implements \link{parse_dts_out_X} for the SIS model
}
