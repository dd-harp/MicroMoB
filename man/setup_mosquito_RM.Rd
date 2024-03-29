% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mosquito_RM.R
\name{setup_mosquito_RM}
\alias{setup_mosquito_RM}
\title{Setup generalized Ross-Macdonald mosquito model}
\usage{
setup_mosquito_RM(
  model,
  stochastic,
  f = 0.3,
  q = 0.9,
  eip,
  p,
  psi,
  nu = 25,
  M,
  Y,
  Z,
  N = NULL
)
}
\arguments{
\item{model}{an object from \link{make_MicroMoB}}

\item{stochastic}{should the model update deterministically or stochastically?}

\item{f}{the blood feeding rate}

\item{q}{the human blood feeding fraction}

\item{eip}{the Extrinsic Incubation Period (may be time varying see \link{time_varying_parameter})}

\item{p}{daily survival probability (may be time and patch varying see \link{time_patch_varying_parameter})}

\item{psi}{a mosquito dispersal matrix (rows must sum to 1)}

\item{nu}{number of eggs laid per oviposition}

\item{M}{total mosquito density per patch (vector of length \code{p})}

\item{Y}{density of incubating mosquitoes per patch (vector of length \code{p})}

\item{Z}{density of infectious mosquitoes per patch (vector of length \code{p})}

\item{N}{\code{l} by \code{p} matrix describing how eggs from mosquitoes in patches are
distributed amongst aquatic habitats. If \code{NULL} it is the identity matrix of dimension
\code{l}.}
}
\value{
no return value
}
\description{
This is a generalized RM model which allows for time varying EIP and
survival probability. It complies with the mosquito component interface, and may
be simulated deterministically or stochastically.
}
