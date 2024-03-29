% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mosquito_RM.R
\name{compute_Z.RM}
\alias{compute_Z.RM}
\title{Compute density of infective mosquitoes for RM model (\eqn{Z})}
\usage{
\method{compute_Z}{RM}(model)
}
\arguments{
\item{model}{an object from \link{make_MicroMoB}}
}
\value{
a vector of length \code{p} giving the density of infected and infectious mosquitoes in each patch
}
\description{
This method returns \code{Z}.
}
