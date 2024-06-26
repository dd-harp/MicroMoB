% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adult-RM.R
\name{make_parameters_MYZ_RM}
\alias{make_parameters_MYZ_RM}
\title{Make parameters for RM adult mosquito model}
\usage{
make_parameters_MYZ_RM(
  pars,
  EIPname,
  p,
  sigma,
  f,
  q,
  nu,
  eggsPerBatch,
  eip,
  calK
)
}
\arguments{
\item{pars}{a \link{list}}

\item{EIPname}{a string: the class name for the EIP model}

\item{p}{daily mosquito survival}

\item{sigma}{emigration rate}

\item{f}{feeding rate}

\item{q}{human blood fraction}

\item{nu}{oviposition rate, per mosquito}

\item{eggsPerBatch}{eggs laid per oviposition}

\item{eip}{maximum length for extrinsic incubation period}

\item{calK}{mosquito dispersal matrix of dimensions \code{nPatches} by \code{nPatches}}
}
\value{
a \link{list}
}
\description{
Make parameters for RM adult mosquito model
}
