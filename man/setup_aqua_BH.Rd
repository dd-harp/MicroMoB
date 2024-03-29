% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aquatic_BH.R
\name{setup_aqua_BH}
\alias{setup_aqua_BH}
\title{Setup aquatic (immature) mosquito model with Beverton-Holt dynamics}
\usage{
setup_aqua_BH(model, stochastic, molt, surv, K, L)
}
\arguments{
\item{model}{an object from \link{make_MicroMoB}}

\item{stochastic}{should the model update deterministically or stochastically?}

\item{molt}{proportion of immature stages which will mature and emerge as adults each day (may be time and patch varying see \link{time_patch_varying_parameter})}

\item{surv}{daily survival probability (may be time and patch varying see \link{time_patch_varying_parameter})}

\item{K}{carrying capacity (may be time and patch varying see \link{time_patch_varying_parameter})}

\item{L}{initial number of immature mosquitoes}
}
\value{
no return value
}
\description{
A single compartment for all aquatic stages is modeled which
suffers density dependent mortality like the Beverton-Holt model.
}
\details{
All parameters can be passed either as a vector of length equal to \code{l}, a matrix with \code{l} rows
and \code{tmax} columns, or a matrix with \code{l} rows and \code{365} columns.
}
